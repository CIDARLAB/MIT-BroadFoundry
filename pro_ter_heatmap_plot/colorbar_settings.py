import matplotlib
import matplotlib.cm as cm

def colorbar_legend_settings(n_logs, scale):
    # This is a placeholder for more complex custom colormaps.
    # When using a standard colormap, enter it in the heatmap_plot
    color_map = 'jetwz'

    settingsD = {}
    
    # Define a new colormap so that '0' values are white or use a standard colormap
    jetwz = [('white')] + [(cm.jet(i)) for i in xrange(1,256)]

    jetwz_map = matplotlib.colors.LinearSegmentedColormap.from_list('new_map', jetwz, N=256)

    cmapD = {'jetwz': jetwz_map}

    settingsD['color_map'] = cmapD[color_map]

    # Generate scale and tick marks for color code bar
    stD = {4: [[r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$'],[0,1,2,3,4]],
           5: [[r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$'],[0,1,2,3,4,5]]}

    settingsD['scale'] = stD[n_logs][0]
    settingsD['ticks'] = stD[n_logs][1]
    
    # Set proper minor tick marks for scale being used
    minor_tickD = {'log':[0.301029996,0.477121255,0.602059991,0.698970004,0.77815125,0.84509804,0.903089987,0.954242509,
                          1.301029996,1.477121255,1.602059991,1.698970004,1.77815125,1.84509804,1.903089987,1.954242509,
                          2.301029996,2.477121255,2.602059991,2.698970004,2.77815125,2.84509804,2.903089987,2.954242509,
                          3.301029996,3.477121255,3.602059991,3.698970004,3.77815125,3.84509804,3.903089987,3.954242509],
                          #4.301029996,4.477121255,4.602059991,4.698970004],
                   'log_on_linear' :[2,3,4,5,6,7,8,9,1,
                                     20,30,40,50,60,70,80,90,
                                     200,300,400,500,600,700,800,900,
                                     2000,3000,4000,5000,6000,7000,8000,9000,
                                     20000,30000,40000,50000,60000,70000,80000,90000]}
    
    settingsD['minorticks'] = minor_tickD[scale]

    return settingsD
