import numpy

def pro_ter_data_selection(mode, parts, dataD, sort):

    chartD = {}

    selected_data = []

    if mode is 'Promoters':
        chartD['legend_title'] = 'Promoters'
        chartD['xaxis_title'] = 'Terminators'
        
        promoters = parts
        
        terminators = dataD['ylabels']

        for x, pro in enumerate(promoters) :
            tervals = []
            for y, ter in enumerate(terminators) :
                tervals.append(dataD[ter][pro])

            selected_data.append(tervals)

    if mode is 'Terminators':
        chartD['legend_title'] = 'Promoters'
        chartD['xaxis_title'] = 'Terminators'
        
        promoters = dataD['xlabels']
  
        terminators = parts

        for y, ter in enumerate(terminators) :
            provals = []
            for x, pro in enumerate(promoters) :
                provals.append(dataD[ter][pro])

            selected_data.append(provals)

    if sort[0] is 'Y' :

        key_list = numpy.array(selected_data[parts.index(sort[1])])
        sorted_data = []

        if mode is 'Promoters':
            ter_array = numpy.array(terminators)
            inds = key_list.argsort()
            new_terminators = ter_array[inds]
            new_promoters = promoters

            for x, pro in enumerate(new_promoters):
                tervals = []
                for y, ter in enumerate(new_terminators) :
                    tervals.append(dataD[ter][pro])

                sorted_data.append(tervals)

        if mode is 'Terminators':
            pro_array = numpy.array(promoters)
            inds = key_list.argsort()
            new_promoters = pro_array[inds]
            new_terminators = terminators

            for y, ter in enumerate(new_terminators):
                provals = []
                for y, ter in enumerate(new_promoters) :
                    provals.append(dataD[ter][pro])
            
                sorted_data.append(provals)

        data = [sorted_data, new_promoters, new_terminators]
        
    else :
        data = [selected_data, promoters, terminators]

    chartD['data'] = data[0]
    chartD['legend_labels'] = data[1]
    chartD['xaxis_labels'] = data[2]

    return chartD
