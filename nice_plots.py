from matplotlib import rc

def format_plots(usetex=False):

    rc('text', usetex=usetex)
    rc('font', family='serif')
    rc('font', size=12.0)
    rc('axes', linewidth=0.4)
    rc('lines', linewidth=0.4)

    return
