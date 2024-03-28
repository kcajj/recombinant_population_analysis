import numpy as np

def viterbi_algorithm(y, A, B, pi):
    """
        viterbi algorithm
        :param y: observation sequence
        :param A: the transition matrix
        :param B: the emission matrix
        :param pi: the initial probability distribution
    """
    print(y)
    print(A)
    print(B)
    print(pi)

    N = B.shape[0]
    print("N: ",N)
    x_seq = np.zeros([N, 0])
    print("x_seq: ",x_seq)
    print("first observation: ",y[0],B[:, y[0]])
    V = np.log(B[:, y[0]]) + np.log(pi)
    print("V: ",V)

    # forward to compute the optimal value function V
    print("loop start")
    for y_ in y[1:]:
        print("y_: ",y_)
        print("B",B[:, y_])
        print("B tile",np.tile(B[:, y_], reps=[N, 1]).T)
        print("A",A.T)
        print("V",np.tile(V, reps=[N, 1]))#probabilities of the past step
        #probability of past step * transition matrix * emission matrix
        '''
        v=[[in the past step we are in A, in the past step we are in B],
            [in the past step we are in A, in the past step we are in B]]
        A=[[transition from A to A, transition from B to A],
            [transition from A to B, transition from B to B]]
        B=[[emissioin from A, emission from A],
            [emission from B, emission from B]
        _v=[[past step A * stay in A * emission from A, past step B * transition to A * emission from A],
            [past step A * transition to B * emission from B, past step B * stay in B * emission from B]]
        '''
        _V = np.log(np.tile(B[:, y_], reps=[N, 1]).T) + np.log(A) + np.tile(V, reps=[N, 1])
        print("_V: ",_V)
        x_ind = np.argmax(_V, axis=1)#mandatorily chose a maximum for both cases, the one in which we end up in A and the one in which we end up in B
        print("x_ind: ",x_ind)
        print("c_", np.c_[x_ind])
        x_seq = np.hstack([x_seq, np.c_[x_ind]])
        print("x_seq: ",x_seq)
        V = _V[np.arange(N), x_ind]
        print("V: ",V)
        print("")
        print("")
        print("")
        
    print("loop end")
    x_T = np.argmax(V)
    print("x_T: ",x_T)

    # backward to fetch optimal sequence
    print("backward")
    x_seq_opt, i = np.zeros(x_seq.shape[1]+1), x_seq.shape[1]
    print("x_seq_opt: ",x_seq_opt,i)
    prev_ind = x_T
    print(prev_ind)
    while i >= 0:
        print("i: ",i)
        x_seq_opt[i] = prev_ind
        print("x_seq_opt[i]=prev_ind: ",x_seq_opt)
        i -= 1
        print("prev_ind: ",prev_ind)
        print("i: ",i)
        print("x_seq",x_seq)
        prev_ind = x_seq[int(prev_ind), i]
        print("prev_ind",prev_ind)
        print("")
        print("")
        print("")
    return x_seq_opt, np.max(V)

i_p=np.array([0.5, 0.5])

t_p=np.array([[0.9,0.1],
              [0.1,0.9]])
e_p=np.array([[0.7,0.2,0.1],
              [0.7,0.1,0.2]])

seq=np.array([1,1,1,1,2,2,2,2])

prediction, likelihood=viterbi_algorithm(seq,t_p,e_p,i_p)

print(prediction)
print(likelihood)

print("")
print("")
print("test with multiple transition probabilities")

def viterbi_algorithm_silent(y, A, B, pi):
    """
        viterbi algorithm
        :param y: observation sequence
        :param A: the transition matrix
        :param B: the emission matrix
        :param pi: the initial probability distribution
    """
    N = B.shape[0]
    x_seq = np.zeros([N, 0])
    V = np.log(B[:, y[0]]) + np.log(pi)
    
    # forward to compute the optimal value function V
    for y_ in y[1:]:
        _V = np.log(np.tile(B[:, y_], reps=[N, 1]).T) + np.log(A.T) + np.tile(V, reps=[N, 1])
        x_ind = np.argmax(_V, axis=1)
        x_seq = np.hstack([x_seq, np.c_[x_ind]])
        V = _V[np.arange(N), x_ind]
    x_T = np.argmax(V)

    # backward to fetch optimal sequence
    x_seq_opt, i = np.zeros(x_seq.shape[1]+1), x_seq.shape[1]
    prev_ind = x_T
    while i >= 0:
        x_seq_opt[i] = prev_ind
        i -= 1
        prev_ind = x_seq[int(prev_ind), i]
    return x_seq_opt, np.max(V)

seq=np.array([1,1,1,1,1,1,1,1,1,0,0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1])
#probability of transition:
probabilities={0.001:0,0.01:0,0.1:0,0.5:0}

for p in probabilities:
    t_p=np.array([[1-p,p],
                  [p,1-p]])
    prediction, likelihood=viterbi_algorithm_silent(seq,t_p,e_p,i_p)
    print(prediction)
    print(likelihood)
    probabilities[p]=likelihood
    print("")

print(probabilities)
print(max(probabilities,key=probabilities.get))

from matplotlib import pyplot as plt

plt.plot(list(probabilities.keys()),list(probabilities.values()))
plt.show()