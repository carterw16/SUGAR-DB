import pickle

def load_pickle(file_path):
    with open(file_path, 'rb') as file:
        return pickle.load(file)


if __name__ == "__main__":
    data = load_pickle('multi_outputs.pkl')
    print("Data loaded from pickle:", data['P_g']['slack'])