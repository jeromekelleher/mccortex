"""
Small script for development.
"""
import _cortex

def main():
    print("main")
    while True: 
        g = _cortex.Graph("A_gen.ctx")
    print(g.get_kmer_size())
    print(g.get_num_cols())
    print(g.get_num_edge_cols())
    print(g.get_capacity())

if __name__ == "__main__":
    main()
