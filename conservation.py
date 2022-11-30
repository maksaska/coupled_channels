from fractions import Fraction


def projection(J):
    result = []

    for i in range(int(2*J) + 1):
        result.append(-J + i)

    return result


def relevant_particle(name: str):
    result = {}

    result['check'] = False
    result['charge'] = 10
    result['strange'] = 10
    result['J'] = Fraction()
    result['I'] = Fraction()
    result['I3'] = Fraction()
    result['mass'] = 0.0
    result['parity'] = +2

    if name == 'pi+':
        result['check'] = True
        result['charge'] = 1
        result['strange'] = 0
        result['J'] = Fraction(0, 1)
        result['I'] = Fraction(1, 1)
        result['I3'] = Fraction(1, 1)
        result['mass'] = 0.13957
        result['parity'] = -1
    if name == 'pi0':
        result['check'] = True
        result['charge'] = 0
        result['strange'] = 0
        result['J'] = Fraction(0, 1)
        result['I'] = Fraction(1, 1)
        result['I3'] = Fraction(0, 1)
        result['mass'] = 0.13498
        result['parity'] = -1
    if name == 'pi-':
        result['check'] = True
        result['charge'] = -1
        result['strange'] = 0
        result['J'] = Fraction(0, 1)
        result['I'] = Fraction(1, 1)
        result['I3'] = Fraction(-1, 1)
        result['mass'] = 0.13957
        result['parity'] = -1
    if name == 'p':
        result['check'] = True
        result['charge'] = 1
        result['strange'] = 0
        result['J'] = Fraction(1, 2)
        result['I'] = Fraction(1, 2)
        result['I3'] = Fraction(1, 2)
        result['mass'] = 0.93827
        result['parity'] = +1
    if name == 'n':
        result['check'] = True
        result['charge'] = 0
        result['strange'] = 0
        result['J'] = Fraction(1, 2)
        result['I'] = Fraction(1, 2)
        result['I3'] = Fraction(-1, 2)
        result['mass'] = 0.93957
        result['parity'] = +1
    if name == 'eta':
        result['check'] = True
        result['charge'] = 0
        result['strange'] = 0
        result['J'] = Fraction(0, 1)
        result['I'] = Fraction(0, 1)
        result['I3'] = Fraction(0, 1)
        result['mass'] = 0.548
        result['parity'] = -1
    if name == 'K+':
        result['check'] = True
        result['charge'] = 1
        result['strange'] = 1
        result['J'] = Fraction(0, 1)
        result['I'] = Fraction(1, 2)
        result['I3'] = Fraction(1, 2)
        result['mass'] = 0.494
        result['parity'] = -1
    if name == 'K0':
        result['check'] = True
        result['charge'] = 0
        result['strange'] = 1
        result['J'] = Fraction(0, 1)
        result['I'] = Fraction(1, 2)
        result['I3'] = Fraction(-1, 2)
        result['mass'] = 0.498
        result['parity'] = -1
    if name == 'K-':
        result['check'] = True
        result['charge'] = -1
        result['strange'] = 1
        result['J'] = Fraction(0, 1)
        result['I'] = Fraction(1, 2)
        result['I3'] = Fraction(-1, 2)
        result['mass'] = 0.494
        result['parity'] = -1
    if name == 'L':
        result['check'] = True
        result['charge'] = 0
        result['strange'] = -1
        result['J'] = Fraction(1, 2)
        result['I'] = Fraction(0, 1)
        result['I3'] = Fraction(0, 1)
        result['mass'] = 1.1157
        result['parity'] = +1
    if name == 'S+':
        result['check'] = True
        result['charge'] = +1
        result['strange'] = -1
        result['J'] = Fraction(1, 2)
        result['I'] = Fraction(1, 1)
        result['I3'] = Fraction(1, 1)
        result['mass'] = 1.1894
        result['parity'] = +1
    if name == 'S0':
        result['check'] = True
        result['charge'] = 0
        result['strange'] = -1
        result['J'] = Fraction(1, 2)
        result['I'] = Fraction(1, 1)
        result['I3'] = Fraction(0, 1)
        result['mass'] = 1.1926
        result['parity'] = +1
    if name == 'S-':
        result['check'] = True
        result['charge'] = -1
        result['strange'] = -1
        result['J'] = Fraction(1, 2)
        result['I'] = Fraction(1, 1)
        result['I3'] = Fraction(-1, 1)
        result['mass'] = 1.1975
        result['parity'] = +1

    return result


def channel_check(name_init: str, name_fin: str):
    result = {}

    result['B_init'] = 'None'
    result['M_init'] = 'None'
    result['B_fin'] = 'None'
    result['M_fin'] = 'None'
    result['path_init'] = 'None'
    result['path_fin'] = 'None'

    B_list = ['p', 'n', 'S+', 'S0', 'S-', 'L']
    M_list = ['pi+', 'pi0', 'pi-', 'eta', 'K+', 'K0', 'K-']

    if name_init in ['pi+p', 'pi0p', 'pi-p', 'pi+n', 'pi0n', 'pi-n']:
        result['path_init'] = 'piN'
    elif name_init in ['etap', 'etan']:
        result['path_init'] = 'etaN'
    elif name_init in ['K+S+', 'K0S+', 'K-S+', 'K+S0', 'K0S0', 'K-S0', 'K+S-', 'K0S-', 'K-S-']:
        result['path_init'] = 'KS'
    elif name_init in ['K+L', 'K0L', 'K-L']:
        result['path_init'] = 'KL'

    if name_fin in ['pi+p', 'pi0p', 'pi-p', 'pi+n', 'pi0n', 'pi-n']:
        result['path_fin'] = 'piN'
    elif name_fin in ['etap', 'etan']:
        result['path_fin'] = 'etaN'
    elif name_fin in ['K+S+', 'K0S+', 'K-S+', 'K+S0', 'K0S0', 'K-S0', 'K+S-', 'K0S-', 'K-S-']:
        result['path_fin'] = 'KS'
    elif name_fin in ['K+L', 'K0L', 'K-L']:
        result['path_fin'] = 'KL'

    B_list = ['p', 'n', 'S+', 'S0', 'S-', 'L']
    M_list = ['pi+', 'pi0', 'pi-', 'eta', 'K+', 'K0', 'K-']

    for meson in M_list:
        if meson in name_init:
            result['M_init'] = meson
            for barion in B_list:
                if barion in name_init.replace(meson, ''):
                    result['B_init'] = barion
        if meson in name_fin:
            result['M_fin'] = meson
            for barion in B_list:
                if barion in name_fin.replace(meson, ''):
                    result['B_fin'] = barion

    return result
