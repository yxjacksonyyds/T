from test import test
import multiprocessing
from multiprocessing import Process, Queue
from contextlib import contextmanager
import sys

@contextmanager
def timeout_managed(func, timeout, *args, **kwargs):
    queue = Queue()
    def wrapper():
        try:
            result = func(*args, **kwargs)
            queue.put(result)
        except Exception as e:
            queue.put(e)
    
    p = Process(target=wrapper)
    p.start()
    p.join(timeout)
    if p.is_alive():
        p.terminate()
        p.join()
        yield None  # Timeout returns None
    else:
        result = queue.get()
        if isinstance(result, Exception):
            raise result
        else:
            yield result

from rdkit import Chem

def is_valid_smiles(smiles):
    """
    Determine if the given SMILES string is valid.

    :param smiles: str, SMILES string
    :return: bool, True if the SMILES string is valid, False otherwise
    """
    try:
        # Try to create a molecule object from the SMILES string
        molecule = Chem.MolFromSmiles(smiles)
        # If molecule creation succeeds and is not None, the SMILES is valid
        return molecule is not None
    except:
        # If any exception occurs during creation, SMILES is invalid
        return False

def get_score(smiles, pref_name):
    t = test()
    qed, sa, affinity = t.cauculate(smiles, pref_name)
    return qed, sa, affinity

def avg(float_list):
    total_sum = sum(float_list)
    count = len(float_list)
    average = total_sum / count
    return average

def deal_file(path):
    Qed = []
    Affinity = []
    SA = []
    time = 0
    valid = 0
    with open(path, 'r') as file:
        for line in file:
            try:
                print("Processing line...")
                sentence = line.split(' || ')
                smiles = sentence[0]
                pref_name = sentence[1][:-1]
                
                if is_valid_smiles(smiles):
                    with timeout_managed(get_score, 120, smiles, pref_name) as result:
                        if result is None:
                            print(f"Timeout for {smiles}, skipping...")
                            continue
                        qed, sa, affinity = result
                        if affinity < -5:
                            Qed.append(qed)
                            SA.append(sa)
                            Affinity.append(affinity)
                            valid = valid + 1
                            print(time, avg(Qed), avg(Affinity), avg(SA), valid)
                else:
                    print('Invalid SMILES')
            except Exception as e:
                print(f"Error processing line: {str(e)}")
                continue
            
            time = time + 1
    return avg(Qed), avg(Affinity), avg(SA)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_file_path>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    print(f"Processing file: {input_file}")
    result = deal_file(input_file)
    print("Final results:")
    print(f"Average QED: {result[0]}")
    print(f"Average Affinity: {result[1]}")
    print(f"Average SA: {result[2]}")
