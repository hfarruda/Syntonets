import numpy as np
import syntonets

#Defining scales
scale_type = 'equal_temperament'

interval2note = [ "C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]
n_repetitionsts = 9
size = len(interval2note)
interval2note = [interval2note[i%size]+str(int(i/size)+1) for i in range(n_repetitionsts*size)]
intervals = [i for i in range(n_repetitionsts*size)]

ratio_to_fundamental_just_scale = [1.,
                               25./24.,
                               9./8.,
                               6./5.,
                               5./4.,
                               4./3.,
                               45./32.,
                               3./2.,
                               8./5.,
                               5./3.,
                               9./5.,
                               15./8.]

ratio_to_fundamental_pythagorean = [1.,
                                 256./243.,
                                 9./8.,
                                 32./27.,
                                 81./64.,
                                 4./3.,
                                 729./512.,
                                 3./2.,
                                 128./81.,
                                 27./16.,
                                 16./9.,
                                 243./128.]

ratio_to_fundamental_meantone = [1.,
                              1.0449,
                              1.1180,
                              1.1963,
                              1.2500,
                              1.3375,
                              1.3975,
                              1.4953,
                              1.5625,
                              1.6719,
                              1.7889,
                              1.8692 ]   

ratio_to_fundamental_werckmeister = [1.,
                                  256./243.,
                                  64./81. * np.sqrt(2),
                                  32./27.,
                                  256./243. * np.power(2,1/4),
                                  4./3.,
                                  1024./729.,
                                  8./9.  * np.power(8,1/4),
                                  128./81.,
                                  1024./729. * np.power(2,1/4),
                                  16./9.,
                                  128./81.  * np.power(2,1/4)]




c = np.power(2.,1./12.)
f1 = (440./np.power(c,9))/8.
print ("F1:", f1)

interval2frequency = []
if scale_type == 'just':
    scale = []
    for i in range(n_repetitionsts):
        scale += (np.array(ratio_to_fundamental_just_scale)*float(i+1)).tolist()
    interval2frequency = [f1 * scale[i] for i in intervals]
elif scale_type == 'pythagorean':
    scale = []
    for i in range(n_repetitionsts):
        scale += (np.array(ratio_to_fundamental_pythagorean)*float(i+1)).tolist()
    interval2frequency = [f1 * scale[i] for i in intervals]
elif scale_type == 'meantone':
    scale = []
    for i in range(n_repetitionsts):
        scale += (np.array(ratio_to_fundamental_meantone)*float(i+1)).tolist()
    interval2frequency = [f1 * scale[i] for i in intervals]
elif scale_type == 'werckmeister':
    scale = []
    for i in range(n_repetitionsts):
        scale += (np.array(ratio_to_fundamental_werckmeister)*float(i+1)).tolist()
    interval2frequency = [f1 * scale[i] for i in intervals]
elif scale_type == 'equal_temperament':
    interval2frequency = [f1 * np.power(c,i) for i in intervals]


g = syntonets.create_network(interval2note, interval2frequency, beta = 1, alpha = 0.2, number_of_edges = 500, giant = True, with_colors = True)
syntonets.visualize(g, './test')
