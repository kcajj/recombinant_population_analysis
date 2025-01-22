import numpy as np
from viterbi import viterbi_algorithm
import matplotlib.pyplot as plt

'''
generate sequences on the basis of given probabilties
'''

initial_p={"A":0.5,"B":0.5}

transition_p_froma={"A":0.99996,"B":0.00004}
transition_p_fromb={"A":0.00004,"B":0.99996}

emission_p_froma={".":0.972,"a":0.027,"b":0.001}
emission_p_fromb={".":0.972,"a":0.001,"b":0.027}

ip_np=np.array(list(initial_p.values()))
tp_np=np.array([list(transition_p_froma.values()),
                list(transition_p_fromb.values())])
ep_np=np.array([list(emission_p_froma.values()),
                list(emission_p_fromb.values())])

accuracies=[]

for n in range(10000,150000,10000):

    for i in range(10):
        
        hid_seq=np.random.choice(list(initial_p.keys()), p=list(initial_p.values()))

        for i in range(n):
            if hid_seq[-1]=="A":
                step=np.random.choice(list(transition_p_froma.keys()), p=list(transition_p_froma.values()))
                hid_seq+=step
            else:
                step=np.random.choice(list(transition_p_fromb.keys()), p=list(transition_p_fromb.values()))
                hid_seq+=step

        emission_seq=""

        for hidden_state in hid_seq:
            if hidden_state=="A":
                emission=np.random.choice(list(emission_p_froma.keys()), p=list(emission_p_froma.values()))
                emission_seq+=emission
            else:
                emission=np.random.choice(list(emission_p_fromb.keys()), p=list(emission_p_fromb.values()))
                emission_seq+=emission

        hidden_index=[]
        for h in hid_seq:
            if h=="A":
                hidden_index.append(0)
            else:
                hidden_index.append(1)
        emission_index=[]
        for e in emission_seq:
            if e==".":
                emission_index.append(0)
            elif e=="a":
                emission_index.append(1)
            else:
                emission_index.append(2)

        prediction, log_lik=viterbi_algorithm(emission_index, tp_np, ep_np, ip_np)

        predicted_hidden_index = prediction.tolist()
        predicted_hidden_index = [int(x) for x in predicted_hidden_index]

        accuracy = sum([1 for i in range(len(hidden_index)) if hidden_index[i] == predicted_hidden_index[i]]) / len(hidden_index)
        accuracies.append(accuracy)

        print(n,accuracy)

print(accuracies)

print(np.mean(accuracies))

'''
[0.9976002399760024, 1.0, 1.0, 1.0, 1.0, 0.9948005199480052, 1.0, 1.0, 1.0, 1.0, 0.9958002099895005, 1.0, 1.0, 0.9935503224838758, 1.0, 0.9957502124893756, 1.0, 1.0, 0.9945502724863757, 0.9943502824858758, 0.9981667277757408, 0.9985667144428519, 0.999600013332889, 0.9998333388887037, 0.9960667977734076, 0.9988333722209259, 0.9991333622212593, 0.9993000233325556, 0.9993333555548148, 1.0, 1.0, 0.9978750531236719, 1.0, 1.0, 0.9984750381240469, 0.9986750331241719, 0.995800104997375, 0.999800004999875, 1.0, 1.0, 1.0, 0.999220015599688, 0.998640027199456, 0.999840003199936, 1.0, 0.996820063598728, 0.99970000599988, 1.0, 1.0, 1.0, 1.0, 0.9991833469442176, 1.0, 0.9980833652772454, 0.9998833352777454, 0.9962000633322778, 0.9999333344444259, 0.9986666888885185, 0.9980833652772454, 0.998683355277412, 1.0, 0.9993714375508921, 1.0, 0.9970857559177726, 0.9994000085713062, 0.9988857302038542, 0.9991857259182012, 0.9975714632648105, 0.9993428665304781, 0.9948857873458951, 1.0, 0.9990875114061074, 0.9999125010937363, 0.9987875151560606, 0.9991125110936113, 0.9941750728115899, 0.9986375170310371, 0.9998000024999687, 0.9983000212497344, 0.9989250134373321, 1.0, 0.9959444895056722, 0.9988000133331851, 0.9980444661725981, 0.9983889067899245, 0.9995444495061167, 0.9968444795057833, 0.9995444495061167, 0.9995333385184609, 0.9960889323451961, 0.998720012799872, 0.998460015399846, 0.997130028699713, 0.999650003499965, 0.999150008499915, 0.997930020699793, 0.998430015699843, 0.997260027399726, 0.997030029699703, 0.996380036199638, 0.9990454632230616, 0.999709093553695, 0.996336396941846, 0.9995363678512014, 0.9987818292560977, 0.9989454641321442, 0.9993000063635785, 0.9990454632230616, 0.9986454668593922, 0.9984091053717693, 0.9996833359722003, 0.9981250156248698, 0.9984750127082275, 0.9970083582636812, 0.9994166715277373, 0.9995083374305214, 0.9970333580553495, 0.9997916684027633, 0.9981416821526488, 0.9989250089582586, 0.999346158875701, 0.9985230882839363, 0.9982846285797802, 0.9979538618933701, 0.9985461650294998, 0.9988307782247828, 0.9995230805916877, 0.999630772070984, 0.9984846270413305, 0.9961231067453328, 0.9975285890815065, 0.9984500110713495, 0.9978357297447875, 0.9961143134691895, 0.9977643016835595, 0.9971071635202605, 0.9983500117856301, 0.9978214441325419, 0.9960643138263298, 0.9990142927550517]
0.9985664119945016
'''