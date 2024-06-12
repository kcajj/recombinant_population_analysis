import numpy as np
from viterbi import viterbi_algorithm
import matplotlib.pyplot as plt

'''
generate sequences on the basis of given probabilties
'''

initial_p={"A":0.5,"B":0.5}

transition_p_froma={"A":0.99999,"B":0.00001}
transition_p_fromb={"A":0.00001,"B":0.99999}

emission_p_froma={".":0.967,"a":0.03,"b":0.003}
emission_p_fromb={".":0.967,"a":0.003,"b":0.03}

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
[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9991000899910009, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.991700414979251, 1.0, 1.0, 1.0, 1.0, 0.9957502124893756, 1.0, 1.0, 1.0, 1.0, 0.9999333355554815, 1.0, 1.0, 0.9968001066631113, 1.0, 0.9995333488883704, 1.0, 0.9961750956226094, 0.9997250068748281, 0.997000074998125, 1.0, 0.9984750381240469, 1.0, 0.9981000474988125, 0.9978250543736407, 1.0, 0.99910001799964, 1.0, 1.0, 1.0, 1.0, 0.999720005599888, 0.998940021199576, 1.0, 0.999560008799824, 1.0, 1.0, 0.9998000033332778, 1.0, 0.9997166713888102, 1.0, 1.0, 0.9989833502774954, 1.0, 0.9997666705554907, 1.0, 0.9993285810202711, 1.0, 0.999600005714204, 0.9995428636733761, 1.0, 1.0, 0.999800002857102, 1.0, 1.0, 1.0, 0.999737503281209, 1.0, 1.0, 0.9990625117186035, 1.0, 1.0, 1.0, 0.9999625004687441, 0.9999125010937363, 1.0, 0.9997777802468861, 0.9997888912345418, 1.0, 0.9998888901234431, 1.0, 1.0, 1.0, 0.9993000077776913, 0.9993777846912812, 1.0, 1.0, 0.999610003899961, 0.999580004199958, 0.999970000299997, 0.999840001599984, 0.99940000599994, 1.0, 0.999650003499965, 1.0, 1.0, 0.9999363642148708, 1.0, 1.0, 0.9998090926446124, 1.0, 0.9994909137189661, 0.9998000018181653, 0.9990909173552968, 0.9999636366942118, 0.9994545504131781, 1.0, 0.9997416688194265, 0.9986750110415746, 0.9991666736110533, 0.9985000124998958, 0.9992083399305006, 0.9990583411804902, 1.0, 1.0, 0.9999000008333264, 0.9979000161537219, 0.9978307859170314, 1.0, 0.9997461557988016, 1.0, 0.999861539526619, 1.0, 0.9993615433727433, 0.999300005384574, 0.9987153944969654, 0.9991642916836309, 0.9988857222448411, 1.0, 0.9990357211734202, 0.9999000007142806, 0.9997214305612103, 0.999921429132649, 0.9972357340304712, 1.0, 1.0]
0.9995514979541436
'''