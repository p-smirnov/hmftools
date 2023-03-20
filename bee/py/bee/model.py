import torch
from torch import nn
from torch.utils.data import DataLoader
from torchvision import datasets
from torchvision.transforms import ToTensor

# best model = 0.0129

NUM_LINEAR_INPUTS = 6
NUM_FCS = 5
FC_HIDDEN_NODES = 100
DROPOUT_P = 0.5
NUM_CHANNELS = 3

'''
        self.conv1 = nn.Conv2d(1, 32, 3, 1)
        self.conv2 = nn.Conv2d(32, 64, 3, 1)
        self.dropout1 = nn.Dropout(0.25)
        self.dropout2 = nn.Dropout(0.5)
        self.fc1 = nn.Linear(9216, 128)
        self.fc2 = nn.Linear(128, 10)

in_channels, out_channels, kernel_size, stride=1, padding=0, dilation=1, groups=1, bias=True, padding_mode='zeros', device=None, dtype=None

'''

# Define model
class NeuralNetwork(nn.Module):
    def __init__(self, read_length):
        super().__init__()
        
        KERNEL_SIZE = 3
        CONV_OUT_CHANNELS = 6
        NUM_CONV_LAYERS = 3
        
        self.conv_stack = nn.Sequential(
            nn.BatchNorm1d(NUM_CHANNELS), # TODO: does this even make sense?
            
            nn.Conv1d(in_channels=NUM_CHANNELS, out_channels=CONV_OUT_CHANNELS, kernel_size=KERNEL_SIZE),
            nn.ReLU(),
            nn.BatchNorm1d(CONV_OUT_CHANNELS),
            nn.Dropout(DROPOUT_P),
            
            nn.Conv1d(in_channels=CONV_OUT_CHANNELS, out_channels=CONV_OUT_CHANNELS, kernel_size=KERNEL_SIZE),
            nn.ReLU(),
            nn.BatchNorm1d(CONV_OUT_CHANNELS),
            nn.Dropout(DROPOUT_P),
            
            nn.Conv1d(in_channels=CONV_OUT_CHANNELS, out_channels=CONV_OUT_CHANNELS, kernel_size=KERNEL_SIZE),
            nn.ReLU(),
            nn.BatchNorm1d(CONV_OUT_CHANNELS),
            nn.Dropout(DROPOUT_P),
        )
        
        # initialise out layer weight
        out_layer = nn.Linear(FC_HIDDEN_NODES, read_length)
        
        # set initials bias such sigmoid(bias) ~ 0.001
        out_layer.bias.data.fill_(-5)
        
        # we have gone through 3 conv layers of kernel size 5
        # we need to calculate the size of the input layer
        output_length = read_length
        for i in range(3):
            output_length = __class__.calc_conv1d_output_length(output_length, KERNEL_SIZE)
        fc_in_size = output_length * CONV_OUT_CHANNELS + NUM_LINEAR_INPUTS
        
        # TODO: check that using sequential dropout is impacted by train / eval
        self.fc_stack = nn.Sequential(
            nn.Linear(fc_in_size, FC_HIDDEN_NODES),
            nn.ReLU(),
            nn.BatchNorm1d(FC_HIDDEN_NODES),
            nn.Dropout(DROPOUT_P),
            nn.Linear(FC_HIDDEN_NODES, FC_HIDDEN_NODES),
            nn.ReLU(),
            nn.BatchNorm1d(FC_HIDDEN_NODES),
            nn.Dropout(DROPOUT_P),
            out_layer,
            nn.Sigmoid() # we can consider using BCEWithLogitsLoss and remove this, but seems a little hairy
        )

    def forward(self, seq_x, data_x):
        
        #print(seq_x.shape)
        x = self.conv_stack(seq_x)
        x = torch.flatten(x, 1)
        
        # concat the sequential input with the linear data
        x = torch.cat((x, data_x), dim=1)
        
        # uncomment follow to work out how big the first linear layer needs to be
        # print(x.shape)
        
        x = self.fc_stack(x)
        return x
    
    # input Cin,Lin, output Cout, Lout
    # https://pytorch.org/docs/stable/generated/torch.nn.Conv1d.html#torch.nn.Conv1d
    @staticmethod
    def calc_conv1d_output_length(length_in, kernel_size, stride=1, padding=0, dilation=1):
        return int((length_in + 2 * padding - dilation * (kernel_size - 1) - 1) / stride + 1)
