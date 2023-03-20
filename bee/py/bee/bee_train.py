import modin.pandas as pd
import ray

import numpy as np
import os, sys, time

import torch
from torch import nn
import torch.utils.data as data_utils
from torchvision import datasets
from torchvision.transforms import ToTensor

import matplotlib.pyplot as plt

import logging

logger = logging.getLogger(__name__)

logging.basicConfig(stream=sys.stdout,
                    format='%(asctime)s.%(msecs)03d | %(threadName)s | %(levelname)s | %(name)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)

logger.setLevel(logging.DEBUG)

# Get cpu or gpu device for training.
#device = "cuda" if torch.cuda.is_available() else "mps" if torch.backends.mps.is_available() else "cpu"

# use cpu as it is much faster for this purpose
device = "cpu"
print(f"Using {device} device")


from torch.utils.data import Dataset

# read bases, base qual, match ref (future: maybe insert)
NUM_CHANNELS = 3

# convert a df to tensor to be used in pytorch
def df_to_tensor(df):
    return torch.from_numpy(df.values).float().to(device)

class BeeDataset(data_utils.Dataset):

    def __init__(self, features_df, targets_df, read_length, transform=None):
        # separate sequence and other data
        self.len = len(features_df)
        
        # do some checks to make sure we got the correct columns
        base_columns = [x for x in features_df.columns if x.startswith("base_") ]
        qual_columns = [x for x in features_df.columns if x.startswith("qual_") ]
        ref_base_columns = [x for x in features_df.columns if x.startswith("refBase_") ]
        
        if len(base_columns) != read_length:
            raise RuntimeError(f"number of base column({len(base_columns)}) != read_length({read_length})")

        if len(qual_columns) != read_length:
            raise RuntimeError(f"number of qual column({len(qual_columns)}) != read_length({read_length})")

        if len(ref_base_columns) != read_length:
            raise RuntimeError(f"number of ref base column({len(ref_base_columns)}) != read_length({read_length})")
            
        seq_tensor = df_to_tensor(features_df[base_columns + qual_columns + ref_base_columns])
        
        # reshape the sequential features into multiple channels
        # todo: check if this is done correctly?
        self.seq_inputs = torch.reshape(seq_tensor, (len(features_df), NUM_CHANNELS, read_length))
        self.linear_inputs = df_to_tensor(features_df[['isRead1', 'isMapped', 'isMateMapped', 'gcContent', 'mappability', 'mapQuality']])
        
        if len(targets_df.columns) != read_length:
            raise RuntimeError(f"number of targets column({len(targets_df.columns)}) != read_length({read_length})")        
        
        self.targets = df_to_tensor(targets_df)
        self.baseq = torch.exp(df_to_tensor(features_df[qual_columns]))
        self.bqr = torch.exp(df_to_tensor(features_df[[x for x in features_df.columns if x.startswith("bqr_")]]))
        self.transform = transform
        
    def __len__(self):
        return self.len
    
    # we get the
    # 1. sequential inputs
    # 2. linear inputs
    # 3. targets
    # 4. the base qual logP
    # 5. the BQR based logP
    def __getitem__(self, idx):
        return self.seq_inputs[idx], self.linear_inputs[idx], self.targets[idx], self.baseq[idx], self.bqr[idx]
    

def load_csv_input(csv_path):
    start = time.time()

    df = pd.read_csv(os.path.expanduser(csv_path), sep="\t",
                 dtype={"isRead1": np.int8, "isMapped": np.int8, "isMateMapped": np.int8})

    elapsed_sec = (time.time() - start)
    minute = int(elapsed_sec / 60)
    second = round(elapsed_sec % 60)
    logger.info(f"loading dataset took {minute}m {second}s")
    return df

#df[df.apply(lambda x: "N" in x["readString"], axis=1)]

def load_feather_input(feather_path):
    start = time.time()

    df = pd.read_feather(os.path.expanduser(feather_path))

    elapsed_sec = (time.time() - start)
    minute = int(elapsed_sec / 60)
    second = round(elapsed_sec % 60)
    logger.info(f"loading dataset took {minute}m {second}s")
    return df

# from the input pandas dataframe create the train / test dataloaders requied for
# pytorch training.
# return (train_dataloader, test_dataloader)
def create_dataloader(df, read_length):
    
    # Using Skicit-learn to split data into training and testing sets
    from sklearn.model_selection import train_test_split

    features_df = df.drop(columns=[c for c in df.columns if c.startswith("baseErr_")])

    targets_df = df[[c for c in df.columns if c.startswith("baseErr_")]]
    targets_df
    
    # Split the data into training and testing sets
    train_features, test_features, train_labels, test_labels = train_test_split(features_df, targets_df, test_size = 0.25, random_state = 42)

    logger.info(f"train size: {len(train_features)}, test size: {len(test_features)}")
    
    # batch size of 100 is unstable, so far 500 looks best
    batch_size = 500

    # Create data loaders.
    train = BeeDataset(train_features, train_labels, read_length)
    train_dataloader = data_utils.DataLoader(train, batch_size=batch_size, shuffle=False)

    test = BeeDataset(test_features, test_labels, read_length)
    test_dataloader = data_utils.DataLoader(test, batch_size=batch_size, shuffle=False)

    #for seq_x, linear_x, y, bq, bqr in test_dataloader:
    #    print(f"Shape of seq x [N, C, H]: {seq_x.shape}")
    #    print(f"Shape of linear x [N, C]: {linear_x.shape}")
    #    print(f"Shape of y: {y.shape} {y.dtype}")
    #    break

    return train_dataloader, test_dataloader

    
def train_model(model, df, epochs, read_length):
    train_dataloader, test_dataloader = create_dataloader(df, read_length)
    
    loss_fn = nn.BCELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    stats_list = []
    for t in range(epochs):
        logger.info(f"Epoch {t+1}\n-------------------------------")
        epoch_stats = EpochStats()
        stats_list.append(epoch_stats)
        train(train_dataloader, model, loss_fn, optimizer, epoch_stats)
        test(test_dataloader, model, loss_fn, epoch_stats)
    logger.info("Done!")

    # create a df for the epoch data
    epoch_df = pd.DataFrame()
    epoch_df["epoch"] = [i + 1 for i in range(epochs)]
    epoch_df["trainLoss"] = [ s.train_loss for s in stats_list ]
    epoch_df["testLoss"] = [ s.test_loss for s in stats_list ]
    epoch_df["bqLoss"] = [ s.bq_loss for s in stats_list ]
    epoch_df["bqrLoss"] = [ s.bqr_loss for s in stats_list ]
    epoch_df["trueErr"] = [ s.true_err for s in stats_list ]
    epoch_df["falseErr"] = [ s.false_err for s in stats_list ]

    ax = plt.subplot()
    epoch_df.plot(ax=ax, x="epoch", y="trainLoss")
    epoch_df.plot(ax=ax, x="epoch", y="testLoss")
    epoch_df.plot(ax=ax, x="epoch", y="bqLoss")
    epoch_df.plot(ax=ax, x="epoch", y="bqrLoss")
    #epoch_df.plot(ax=ax, x="epoch", y="trueErr")
    #epoch_df.plot(ax=ax, x="epoch", y="falseErr")
    
class EpochStats:
    def __init__(self):
        pass

def train(dataloader, model, loss_fn, optimizer, epoch_stats):
    size = len(dataloader.dataset)
    model.train()
    train_loss = 0
    for batch, (seq_x, lin_x, y, bq, bqr) in enumerate(dataloader):
        
        # Compute prediction error
        pred = model(seq_x, lin_x)
        loss = loss_fn(pred, y)

        # Backpropagation
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
        train_loss += loss.item()

        if batch % 100 == 0:
            loss = loss.item()
            current = (batch + 1) * y.size(dim=0)
            logger.debug(f"loss: {loss:>7f}  [{current:>5d}/{size:>5d}]")
            
    epoch_stats.train_loss = train_loss / len(dataloader)
    
def calc_error(pred, target):
    num_target_true = target.sum().item()
    num_target_false = (1 - target).sum().item()

    true_error = (target - pred) * target

    # 1 - target gives the mask required to zero out all true target ones
    false_error = (pred - target) * (1 - target)

    av_true_error = true_error.sum().item() / num_target_true if num_target_true > 0 else float('nan')
    av_false_error = false_error.sum().item() / num_target_false if num_target_false > 0 else float('nan')
    
    #logger.debug(f"pred: {torch.log(pred)}, target: {target}, av true err:{av_true_error}, av false err: {av_false_error}")
    
    return av_true_error, av_false_error
    
def test(dataloader, model, loss_fn, epoch_stats):
    size = len(dataloader.dataset)
    num_batches = len(dataloader)
    model.eval()
    test_loss = 0
    bq_loss = 0
    bqr_loss = 0
    true_err, false_err, num_true, num_false = 0, 0, 0, 0
    with torch.no_grad():
        for seq_x, lin_x, y, bq, bqr in dataloader:
            pred = model(seq_x, lin_x)
            #pred = torch.full(y.shape, 1.0)
            test_loss += loss_fn(pred, y).item()
            bq_loss += loss_fn(bq, y).item()
            bqr_loss += loss_fn(bqr, y).item()
            
            t_err, f_err = calc_error(pred, y)
            if pd.notna(t_err):
                true_err += t_err
                num_true += 1
            if pd.notna(f_err):
                false_err += f_err
                num_false += 1
            #logger.debug(f"true error: {t_err}, false error: {f_err}")
    test_loss /= num_batches
    bq_loss /= num_batches
    bqr_loss /= num_batches
    
    true_err /= num_true
    false_err /= num_false
    
    epoch_stats.test_loss = test_loss
    epoch_stats.bq_loss = bq_loss
    epoch_stats.bqr_loss = bqr_loss
    
    epoch_stats.true_err = true_err
    epoch_stats.false_err = false_err
    logger.debug(f"Test Error: \n Accuracy: true error: {true_err:.4f}, false error: {false_err:.4f}, Avg loss: {test_loss:>8f} \n")
    
def main():
    logger.info(f"starting bee trainer")
    parser = argparse.ArgumentParser(description='Bee trainer')
    parser.add_argument('input_csv', help='input csv file')
    parser.add_argument('outdir', help='output directory')
    parser.add_argument('read_length', help='read length')
    args = parser.parse_args()

    tasks = []


if __name__ == "__main__":
    main()

