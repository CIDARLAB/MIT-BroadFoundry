def pro_ter_data_selection(mode, parts, dataD):

    chartD = {}

    data = []

    if mode is 'Promoters':
        chartD['legend_title'] = 'Promoters'
        chartD['xaxis_title'] = 'Terminators'
        
        promoters = parts
        chartD['legend_labels'] = promoters
        
        terminators = dataD['ylabels']

        chartD['xaxis_labels'] = terminators

        for x, pro in enumerate(promoters) :
            tervals = []
            for y, ter in enumerate(terminators) :

                tervals.append(dataD[ter][pro])

            data.append(tervals)


    if mode is 'Terminators':
        chartD['legend_title'] = 'Promoters'
        chartD['xaxis_title'] = 'Terminators'
        
        promoters = dataD['xlabels']
        chartD['xaxis_labels'] = promoters
        
        terminators = parts
        chartD['legend_labels'] = terminators

        for y, ter in enumerate(terminators) :
            provals = []
            for x, pro in enumerate(promoters) :
                provals.append(dataD[ter][pro])

            data.append(provals)

    chartD['data'] = data

    print data


    return chartD
