import matplotlib.pyplot as plt

def set_style():
    #Direct input 
    plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
    #Options
    params = {'text.usetex' : True,
              'font.size' : 20,
              'font.family' : 'lmodern',
              'text.latex.unicode': True,
              }
    plt.rcParams.update(params) 
    
    return 
