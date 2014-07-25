"""
Small script for development.
"""
import _cortex

def memory_test():
     while True: 
        g = _cortex.Graph("A_gen.ctx")
        kmer_size = g.get_kmer_size()
        num_cols = g.get_num_cols()
        for k in _cortex.KmerIterator(g):
            assert g.contains_kmer(k)
            t = g.get_next_nodes(k, -1)
            for j in range(num_cols):
                c = g.contains_kmer(k, j)
        kmer = "G" * 31
        try:
            t = g.get_next_nodes(kmer, 0)
        except:
            pass
        #break

def main():
    g = _cortex.Graph("A_gen.ctx")
    print(g.get_kmer_size())
    print(g.get_num_cols())
    print(g.get_num_edge_cols())
    print(g.get_capacity())
    #print(len(kmers))
    col = 0
    for k in _cortex.KmerIterator(g):
        assert g.contains_kmer(k)
        t = g.get_next_nodes(k)
        #print(k, t)
        print(k, g.contains_kmer(k, col))
    kmer = "G" * 31
    print(g.contains_kmer(kmer))
    try:
        t = g.get_next_nodes(kmer, 0)
    except:
        pass


if __name__ == "__main__":
    #main()
    memory_test()
