# go-nuc-iter 
- a go implementation the same speed as the RUST. 
- to compare against RUST, i coded this also in RUST.
- golang has no unique function, so implemented a additional string compare. 

```
╭─gauavsablok@gauravsablok ~/Desktop/go/go-nuc-iter
╰─$ go run main.go -h
Finding the origin of the kmer

Usage:
  sequence [flags]

Flags:
  -G, --genome file string    genome file to be analyzed (default "genome file")
  -h, --help                  help for sequence
  -I, --illuminafile string   illumina file to be analyzed (default "illumina file")
  -A, --kmerArgs int          origin kmer (default 10)
  -P, --pacbiofile string     pacbioinsert (default "pacbiofile to be analyzed")
exit status 1
``

Gaurav Sablok
