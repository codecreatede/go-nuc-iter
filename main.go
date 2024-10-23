package main

/*

Author Gaurav Sablok
Universitat Potsdam
Date : 2024-10-23

A golang implementation of the same:

 A kmer origin finding faster than the recent implementation of the recent implementation
 Back to sequences: Find the origin of 𝑘-mers DOI: 10.21105/joss.07066.

 I implemented the rust async programming to index the kmer first over a window size and then
 use that to make the set of the kmers, so that you have less search space and using that to
 search the kmer in the file provided

 it only searchers for the unique hashes and their location. to make it even faster, i am also
 implementing a async programming later today.

 it support genome and short and long illumina reads.

--- golang has no string compare function for the unique hashes, so implemented a iter to compare
the hashes, if you want you can store them as BTrees also.


*/

import (
	"bufio"
	"log"
	"os"
	"strings"

	"github.com/spf13/cobra"
)

var (
	pacbiofile   string
	illuminafile string
	genomefile   string
	kmerArgs     int
)

var rootCmd = &cobra.Command{
	Use:  "sequence",
	Long: "Finding the origin of the kmer",
	Run:  sequenceFunc,
}

func main() {
	if err := rootCmd.Execute(); err != nil {
		log.Fatal(err)
	}
	os.Exit(1)
}

func init() {
	rootCmd.Flags().
		StringVarP(&pacbiofile, "pacbiofile", "P", "pacbiofile to be analyzed", "pacbioinsert")
	rootCmd.Flags().
		StringVarP(&illuminafile, "illuminafile", "I", "illumina file", "illumina file to be analyzed")
	rootCmd.Flags().
		StringVarP(&genomefile, "genome file", "G", "genome file", "genome file to be analyzed")
	rootCmd.Flags().
		IntVarP(&kmerArgs, "kmerArgs", "A", 10, "origin kmer")
}

func sequenceFunc(cmd *cobra.Command, args []string) {
	type pacbiofileID struct {
		id string
	}
	type pacbiofileSeq struct {
		seq string
	}

	type illuminaID struct {
		id string
	}

	type illuminaSeq struct {
		seq string
	}

	type genomeID struct {
		id string
	}

	type genomeSeq struct {
		seq string
	}

	pacbioIDConstruct := []pacbiofileID{}
	pacbioSeqConstruct := []pacbiofileSeq{}

	illuminaIDConstruct := []illuminaID{}
	illuminaSeqConstruct := []illuminaSeq{}

	genomeIDConstruct := []genomeID{}
	genomeSeqConstruct := []genomeSeq{}

	// reading and storing the pacbio file

	fpacbio, err := os.Open(pacbiofile)
	if err != nil {
		log.Fatal(err)
	}
	Opacbio := bufio.NewScanner(fpacbio)
	for Opacbio.Scan() {
		line := Opacbio.Text()
		if strings.HasPrefix(string(line), "@") {
			pacbioIDConstruct = append(pacbioIDConstruct, pacbiofileID{
				id: strings.ReplaceAll(strings.Split(string(line), " ")[0], "@", ""),
			})
		}
		if strings.HasPrefix(string(line), "A") || strings.HasPrefix(string(line), "T") ||
			strings.HasPrefix(string(line), "G") ||
			strings.HasPrefix(string(line), "C)") {
			pacbioSeqConstruct = append(pacbioSeqConstruct, pacbiofileSeq{
				seq: string(line),
			})
		}
	}

	// reading and storing the illumina file

	fillumina, err := os.Open(illuminafile)
	if err != nil {
		log.Fatal(err)
	}
	Oillumina := bufio.NewScanner(fillumina)
	for Oillumina.Scan() {
		line := Oillumina.Text()
		if strings.HasPrefix(string(line), "@") {
			illuminaIDConstruct = append(illuminaIDConstruct, illuminaID{
				id: strings.ReplaceAll(strings.Split(string(line), " ")[0], "@", ""),
			})
		}
		if strings.HasPrefix(string(line), "A") || strings.HasPrefix(string(line), "T") ||
			strings.HasPrefix(string(line), "G") ||
			strings.HasPrefix(string(line), "C)") {
			illuminaSeqConstruct = append(illuminaSeqConstruct, illuminaSeq{
				seq: string(line),
			})
		}
	}

	// reading and storing the genome file

	fgenome, err := os.Open(genomefile)
	if err != nil {
		log.Fatal(err)
	}
	Ogenome := bufio.NewScanner(fgenome)
	for Ogenome.Scan() {
		line := Ogenome.Text()
		if strings.HasPrefix(string(line), ">") {
			genomeIDConstruct = append(genomeIDConstruct, genomeID{
				id: strings.ReplaceAll(strings.Split(string(line), " ")[0], ">", ""),
			})
		}
		if strings.HasPrefix(string(line), "A") || strings.HasPrefix(string(line), "T") ||
			strings.HasPrefix(string(line), "G") ||
			strings.HasPrefix(string(line), "C)") {
			genomeSeqConstruct = append(genomeSeqConstruct, genomeSeq{
				seq: string(line),
			})
		}
	}

	// isolating the sequences

	pacbioIsolate := []string{}
	illuminaIsolate := []string{}
	genomeIsolate := []string{}

	for i := range pacbioSeqConstruct {
		pacbioIsolate = append(pacbioIsolate, pacbioSeqConstruct[i].seq)
	}

	for i := range illuminaSeqConstruct {
		illuminaIsolate = append(illuminaIsolate, illuminaSeqConstruct[i].seq)
	}

	for i := range genomeSeqConstruct {
		genomeIsolate = append(genomeIsolate, genomeSeqConstruct[i].seq)
	}

	// making the hashes

	pacbioHashes := []string{}
	illuminaHashes := []string{}
	genomeHashes := []string{}

	for i := 0; i <= len(pacbioIsolate); i++ {
		for j := 0; j <= len(pacbioIsolate[i])-kmerArgs; j++ {
			pacbioHashes = append(pacbioHashes, pacbioIsolate[i][j:j+kmerArgs])
		}
	}

	for i := 0; i <= len(illuminaIsolate); i++ {
		for j := 0; j <= len(pacbioIsolate[i])-kmerArgs; j++ {
			illuminaHashes = append(illuminaHashes, pacbioIsolate[i][j:j+kmerArgs])
		}
	}

	for i := 0; i <= len(genomeIsolate); i++ {
		for j := 0; j <= len(pacbioIsolate[i])-kmerArgs; j++ {
			genomeHashes = append(genomeHashes, pacbioIsolate[i][j:j+kmerArgs])
		}
	}

	// calling the unqiuehashes here

	uniquePacbio, _ := uniqueHash(pacbioIsolate)
	uniqueillumina, _ := uniqueHash(illuminaIsolate)
	uniqueGenome, _ := uniqueHash(genomeIsolate)

	pacbioFile, err := os.Create("pacbiohashes.txt")
	if err != nil {
		log.Fatal(err)
		defer pacbioFile.Close()
	}
	for i := range uniquePacbio {
		pacbioFile.WriteString(uniquePacbio[i] + "\n")
	}

	illuminaFile, err := os.Create("pacbiohashes.txt")
	if err != nil {
		log.Fatal(err)
		defer illuminaFile.Close()
	}
	for i := range uniqueillumina {
		illuminaFile.WriteString(uniqueillumina[i] + "\n")
	}

	genomeFile, err := os.Create("pacbiohashes.txt")
	if err != nil {
		log.Fatal(err)
		defer genomeFile.Close()
	}
	for i := range uniqueGenome {
		genomeFile.WriteString(uniqueGenome[i] + "\n")
	}

	pacbioIndexStart := []int{}
	pacbioIndexEnd := []int{}

	genomeIndexStart := []int{}
	genomeIndexEnd := []int{}

	illuminaIndexStart := []int{}
	illuminaIndexEnd := []int{}

	for i := range uniquePacbio {
		for j := range pacbioIsolate {
			start := strings.Index(pacbioIsolate[j], uniquePacbio[i])
			end := start + len(uniquePacbio[i])
			pacbioIndexStart = append(pacbioIndexStart, start)
			pacbioIndexEnd = append(pacbioIndexEnd, end)
		}
	}

	for i := range uniqueillumina {
		for j := range illuminaIsolate {
			start := strings.Index(illuminaIsolate[j], uniqueillumina[i])
			end := start + len(uniqueillumina[i])
			illuminaIndexStart = append(illuminaIndexStart, start)
			illuminaIndexEnd = append(illuminaIndexEnd, end)
		}
	}

	for i := range uniqueGenome {
		for j := range genomeIsolate {
			start := strings.Index(genomeIsolate[j], uniqueGenome[i])
			end := start + len(uniqueGenome[i])
			genomeIndexStart = append(genomeIndexStart, start)
			genomeIndexEnd = append(genomeIndexEnd, end)
		}
	}

	pacbioWrite, err := os.Create("pacbioKmerOrigin.txt")
	if err != nil {
		log.Fatal(err)
	}

	for i := range pacbioIndexStart {
		pacbioWrite.WriteString(
			string(pacbioIndexStart[i]) + "\t" + string(pacbioIndexEnd[i]) + "\t" + uniquePacbio[i],
		)
	}

	illuminaWrite, err := os.Create("illuminaKmerOrigin.txt")
	if err != nil {
		log.Fatal(err)
	}

	for i := range illuminaIndexStart {
		illuminaWrite.WriteString(
			string(
				illuminaIndexStart[i],
			) + "\t" + string(
				illuminaIndexEnd[i],
			) + "\t" + uniqueillumina[i],
		)
	}

	genomeWrite, err := os.Create("genomeKmerOrigin.txt")
	if err != nil {
		log.Fatal(err)
	}

	for i := range genomeIndexStart {
		genomeWrite.WriteString(
			string(genomeIndexStart[i]) + "\t" + string(genomeIndexEnd[i]) + "\t" + uniqueGenome[i],
		)
	}
}

// golang has no unique implementation, so a additional unqiue function to make the hashes compare to each other.
func uniqueHash(inputvar []string) ([]string, int) {
	captureHash := inputvar
	captureUnique := []string{}
	var captureLength int

	for i := 0; i <= len(captureHash)-1; i++ {
		if captureHash[i] == captureHash[i+1] {
			continue
		} else {
			captureUnique = append(captureUnique, captureHash[i])
		}
	}
	captureLength += len(captureUnique)

	return captureUnique, captureLength
}
