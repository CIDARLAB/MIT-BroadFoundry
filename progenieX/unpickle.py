import pickle

def maine() :
    depickle(raw_input('File:'))
    
def depickle(input_file) :

    with open(input_file, 'r') as p :
        up = pickle.Unpickler(p)
        up_object = up.load()
    print up_object
    return up_object

if __name__ == "__maine__" :
    maine()
