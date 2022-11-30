from interaction import Interaction as I
from conservation import channel_check as ch
import sys

import pandas as pd


def main():
    channel_init = sys.argv[1]
    channel_fin = sys.argv[2]

    model = I(channel_init, channel_fin)   

    df = pd.DataFrame({'W':[1.535], 'theta':[0.0], 'phi':[0.0]})
    #print(model.diff_cs(df))

    #model.check_diff_cs(1.958)
    #model.check_diff_cs(1.443)
    #model.check_diff_cs(1.73)
    #model.check_diff_cs(1.832)
    #model.check_diff_cs(1.999)
    model.check_diff_cs(1.822)

if __name__ == '__main__':
    main()
