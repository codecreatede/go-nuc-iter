# go-nuc-iter 
- a go implementation the same speed as the RUST. 
- to compare against RUST, a RUSt version with the RUST docker is also present.
- golang has no unique function, so implemented a additional string compare. 

```
╭─gauavsablok@gauravsablok ~/Desktop/go/go-nuc-iter ‹main●›
╰─$ go run main.go -h
Finding the origin of the kmer

Usage:
  sequence [command]

Available Commands:
  completion  Generate the autocompletion script for the specified shell
  genome
  help        Help about any command
  illumina
  pacbio

Flags:
  -h, --help   help for sequence

Use "sequence [command] --help" for more information about a command.
exit status 1
╭─gauavsablok@gauravsablok ~/Desktop/go/go-nuc-iter ‹main●›
╰─$ go run main.go pacbio -h                                                                                                                                        1 ↵
Pacbio finding the original kmer

Usage:
  sequence pacbio [flags]

Flags:
  -h, --help                help for pacbio
  -A, --kmerArgs int        origin kmer (default 10)
  -P, --pacbiofile string   pacbioinsert (default "pacbiofile to be analyzed")
exit status 1
╭─gauavsablok@gauravsablok ~/Desktop/go/go-nuc-iter ‹main●›
╰─$ go run main.go illumina -h                                                                                                                                      1 ↵
illumina finding the original kmer

Usage:
  sequence illumina [flags]

Flags:
  -h, --help                  help for illumina
  -i, --illuminafile string   pacbioinsert (default "pacbiofile to be analyzed")
  -A, --kmerArgs int          origin kmer (default 10)
exit status 1
╭─gauavsablok@gauravsablok ~/Desktop/go/go-nuc-iter ‹main●›
╰─$ go run main.go genome -h                                                                                                                                        1 ↵
Genome finding the original kmer

Usage:
  sequence genome [flags]

Flags:
  -G, --genomefile string   pacbioinsert (default "pacbiofile to be analyzed")
  -h, --help                help for genome
  -A, --kmerArgs int        origin kmer (default 10)
exit status 1

```

- it will produce a SAM compatible format which will give the start, end, kmer, origin of the kmer and all the other details. 
- writing a interconvertible which will convert this into a hash Btree. 

```
0       5       ATACT   ATACTTTAAATTTTAGTTACTATTAT
1       6       TACTT   ATACTTTAAATTTTAGTTACTATTAT
2       7       ACTTT   ATACTTTAAATTTTAGTTACTATTAT
3       8       CTTTA   ATACTTTAAATTTTAGTTACTATTAT
4       9       TTTAA   ATACTTTAAATTTTAGTTACTATTAT
```

Gaurav Sablok
