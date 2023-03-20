import modin.pandas as pd
import ray
import numpy as np
import math

def expand_columns(df, read_length):
    expand_read_string(df, read_length)
    expand_base_qual(df, read_length)
    expand_base_err(df, read_length)
    expand_bqr(df, read_length)
    
def expand_read_string(df, read_length):
    # convert the read string using encoding
    data = {}

    # encode base to spread it out over the values
    # we have to reserve the value for 0 when it do not match
    def encode_base(base):
        if base == "A":
            return 0.4
        if base == "C":
            return 0.6
        if base == "T":
            return 0.8
        if base == "G":
            return 1.0
        else:
            return 0

    # create a column for each base
    df.drop(columns = [ c for c in df.columns if c.startswith("base_")], inplace=True)
    df.drop(columns = [ c for c in df.columns if c.startswith("refBase_")], inplace=True)
    for i in range(read_length):
        a = df['readString'].apply(lambda x: encode_base(x[i]) if len(x) > i else encode_base('N'))
        df[f"base_{i}"] = a.astype(np.float16)

        a = df['refGenomeBases'].apply(lambda x: encode_base(x[i]) if len(x) > i else encode_base('N'))
        df[f"refBase_{i}"] = a.astype(np.float16)

def expand_base_qual(df, read_length):
    
    # assert(1 / math.exp(phred_to_logp(30)) == 1000)

    # expand the base qualities
    df.drop(columns = [ c for c in df.columns if c.startswith("qual_")], inplace=True)
    for i in range(read_length):
        df[f"qual_{i}"] = df['baseQualities'].apply(lambda x: adjusted_phred(x[i]) if len(x) > i else float('nan')).astype(np.float16)

def phred_q(code):
    assert(len(code) == 1)
    return ord(code[0]) - 33

# normalise the number, phred highest is 42 (K)
# however in reality 37 (F) is max we can see
# 2 is the smallest in ours it seems
def normalised_phred(code):
    return phred_q(code) / 42.0

# we change the phred to log(P) error
def adjusted_phred(code):
    return phred_to_logp(phred_q(code))
    #return phred_q(code)

# convert phred to log(P) error
# phred is -10log10(P)
# phred = -10log10(P)
# phred / -10.0 = log10(P)
# log(P) = log10(P) / log10(e) = phred / -10.0 / log10(e)
def phred_to_logp(phred):
    return (phred / -10.0) / math.log10(math.e)

def expand_base_err(df, read_length):
    has_error_base = pd.Series([False] * len(df))

    if "hasErrorBase" in df.columns:
        df.drop(columns = "hasErrorBase", inplace=True)
    df.drop(columns = [ c for c in df.columns if c.startswith("baseErr_")], inplace=True)
    # expand the target column
    for i in range(read_length):
        # not sure whether should use 1 or 0 for columns after max length
        target = df['baseError'].apply(lambda x: int(x[i]) if len(x) > i else 1).astype(np.int8)
        df[f"baseErr_{i}"] = target
        has_error_base |= (target == 1)

    df.insert(0, "hasErrorBase", has_error_base)

def expand_bqr(df, read_length):
    # now expand the BQR columns, make them into log(P) error
    df.drop(columns = [ c for c in df.columns if c.startswith("bqr_")], inplace=True)
    bqrQualSplit = df['bqrQualities'].str.split(',')
    for i in range(read_length):
        df[f"bqr_{i}"] = bqrQualSplit.apply(lambda x: phred_to_logp(float(x[i])) if len(x) > i else float('nan')).astype(np.float16)
