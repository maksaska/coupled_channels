from interaction import Interaction as I
from conservation import channel_check as ch
import pandas as pd
import sys

def main():
    channel_init = sys.argv[1]
    channel_fin = sys.argv[2]
    W = sys.argv[3]

    model = I(channel_init, channel_fin)

    model.check_diff_cs(W, reduced=True, mult=[])

    #df = pd.DataFrame({'W':[1.73], 'theta':[141.0], 'phi':[10.0]})
    #print(model.diff_cs(df, reduced=False, mult=[], check=False))

    # ------- Differential Cross Section Plots (theta) -------

    #model.check_diff_cs(1.958)
    #model.check_diff_cs(1.443, reduced=True, mult=['P11'])
    #model.check_diff_cs(1.443, reduced=True, mult=['P13'])
    #model.check_diff_cs(1.443, reduced=True, mult=['P31'])
    #model.check_diff_cs(1.443, reduced=True, mult=['P33'])
    #model.check_diff_cs(1.73, reduced=True, mult=['S11','S31','P11', 'P13', 'P13', 'P33', 'P31', 'D13', 'D15', 'D33', 'D35', 'F15', 'F17', 'F35'])
    #model.check_diff_cs(1.73, reduced=True, mult=['G17'])
    #model.check_diff_cs(1.73, reduced=True, mult=['P31'])
    #model.check_diff_cs(1.73, reduced=True, mult=['P33'])
    #model.check_diff_cs(1.443, reduced=True, mult=[])
    #model.check_diff_cs(1.832)
    #model.check_diff_cs(1.999)
    
    #model.check_diff_cs(1.73)
    #model.check_diff_cs(1.443)
    #model.check_diff_cs(1.958)

    # ------- Differential Cross Section Plots (W) -------

    #model.check_diff_cs_W(140, reduced=True, mult=['G39'])

    # ------- Differential Cross Section Plots (theta) with different multipoles -------

    show_mult = False

    if show_mult:
        sum = []

        sum_sort = sorted(list(model.Data.keys()), key=lambda state : 100*['S', 'P', 'D', 'F', 'G', 'H', 'I'].index(list(state)[0]) + 10*int(list(state)[2]) + int(list(state)[1]))

        for i in sum_sort:
            model.check_diff_cs(1.73, reduced=True, mult=[i])

        for i in sum_sort:
            sum.append(i)

            model.check_diff_cs(1.73, reduced=True, mult=sum)

    # ------- Differential Cross Section Plots (phi) -------

    #model.check_diff_cs_phi(1.443, 50)

    # ------- Multipoles Interpolation Plots -------

    #model.check_interpolation('D13')
    #model.check_interpolation('D15')
    #model.check_interpolation('D33')
    #model.check_interpolation('D35')
    #model.check_interpolation('F15')
    #model.check_interpolation('F17')
    #model.check_interpolation('F35')
    #model.check_interpolation('F37')
    #model.check_interpolation('G17')
    #model.check_interpolation('G19')
    #model.check_interpolation('P13')
    #model.check_interpolation('P11')
    #model.check_interpolation('P31')
    ##model.check_interpolation('P33')
    #model.check_interpolation('S11')
    #model.check_interpolation('H19')
    #model.check_interpolation('H39')    

if __name__ == '__main__':
    main()
