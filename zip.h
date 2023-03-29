#ifndef PARFAIT_ZIP_H
#define PARFAIT_ZIP_H

typedef struct Data_info{
    int data_num;
    char *data_name;
    int *freq_arr;
    int *file_data;
    int total_freq;
    int data_len;
}Data_info; //data about the fasta file

typedef struct Heap_node{
    struct Heap_node *left_node;
    struct Heap_node *right_node;
    char *ch;
}Heap_node; //A node from heap

typedef struct Element{
    int freq;
    Heap_node *pnode;
}Element; //data of heap node

typedef struct Heap{
    Element *heap;
    int size;
}Heap; //heap

static int power5(int power);
static bool isfafile(char *f); //cheak the file is fasta
static char* read_line(FILE *fp, char st); //read one line of the file and return it as a string
static void read_fafile(FILE *fp, FILE *w_fp, int tmp_len); //read the fasta file
static int chr_to_num(char tmp);
static int str_to_num(const char *str, int tmp_len); // ATGCN to 01234
static char* num_to_chr(int idx, int tmp_len); //01234 to ATGAN
static Element mk_huffman_tree(Data_info data_tmp, char *odd_data_tmp, int tmp_len); //making huffman tree
static void heap_insert(Heap *h, Element node); //insert node to heap
static Element heap_delete(Heap *h); //delete node from heap
static int encoding(Element head, Data_info data_tmp, FILE *w_fp, int tmp_len); //encoding file with huffman tree
static void destroy_heap(Heap_node *node);
static void make_huffman_code(Heap_node *node, int len, char *code, char **code_arr, int tmp_len);
static void write_index_file(Data_info data_tmp);
int main_zip(int argc, char *argv[]);

#endif //PARFAIT_ZIP_H
