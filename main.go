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
	"strconv"
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
}

var pacbioCmd = &cobra.Command{
	Use:  "pacbio",
	Long: "Pacbio finding the original kmer",
	Run:  pacbioFunc,
}

var genomeCmd = &cobra.Command{
	Use:  "genome",
	Long: "Genome finding the original kmer",
	Run:  genomeFunc,
}

var illuminaCmd = &cobra.Command{
	Use:  "illumina",
	Long: "illumina finding the original kmer",
	Run:  illuminaFunc,
}

func main() {
	if err := rootCmd.Execute(); err != nil {
		log.Fatal(err)
	}
	os.Exit(1)
}

func init() {
	pacbioCmd.Flags().
		StringVarP(&pacbiofile, "pacbiofile", "P", "pacbiofile to be analyzed", "pacbioinsert")
	pacbioCmd.Flags().
		IntVarP(&kmerArgs, "kmerArgs", "A", 10, "origin kmer")
	illuminaCmd.Flags().
		StringVarP(&illuminafile, "illuminafile", "i", "pacbiofile to be analyzed", "pacbioinsert")
	illuminaCmd.Flags().
		IntVarP(&kmerArgs, "kmerArgs", "A", 10, "origin kmer")
	genomeCmd.Flags().
		StringVarP(&genomefile, "genomefile", "G", "pacbiofile to be analyzed", "pacbioinsert")
	genomeCmd.Flags().
		IntVarP(&kmerArgs, "kmerArgs", "A", 10, "origin kmer")

	rootCmd.AddCommand(pacbioCmd)
	rootCmd.AddCommand(illuminaCmd)
	rootCmd.AddCommand(genomeCmd)
}

func pacbioFunc(cmd *cobra.Command, args []string) {
	type pacbiofileID struct {
		id string
	}
	type pacbiofileSeq struct {
		seq string
	}
	pacbioIDConstruct := []pacbiofileID{}
	pacbioSeqConstruct := []pacbiofileSeq{}

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

	pacbioIsolate := []string{}

	for i := range pacbioSeqConstruct {
		pacbioIsolate = append(pacbioIsolate, pacbioSeqConstruct[i].seq)
	}

	pacbioHashes := []string{}

	for i := 0; i <= len(pacbioIsolate); i++ {
		for j := 0; j <= len(pacbioIsolate[i])-kmerArgs; j++ {
			pacbioHashes = append(pacbioHashes, pacbioIsolate[i][j:j+kmerArgs])
		}
	}

	uniquePacbio, _ := uniqueHash(pacbioIsolate)

	pacbioFile, err := os.Create("pacbiohashes.txt")
	if err != nil {
		log.Fatal(err)
		defer pacbioFile.Close()
	}
	for i := range uniquePacbio {
		pacbioFile.WriteString(uniquePacbio[i] + "\n")
	}

	pacbioIndexStart := []int{}
	pacbioIndexEnd := []int{}

	for i := range uniquePacbio {
		for j := range pacbioIsolate {
			start := strings.Index(pacbioIsolate[j], uniquePacbio[i])
			end := start + len(uniquePacbio[i])
			pacbioIndexStart = append(pacbioIndexStart, start)
			pacbioIndexEnd = append(pacbioIndexEnd, end)
		}
	}

	pacbioWrite, err := os.Create("pacbioKmerOrigin.txt")
	if err != nil {
		log.Fatal(err)
	}

	for i := range pacbioIndexStart {
		start := strconv.Itoa(pacbioIndexStart[i])
		end := strconv.Itoa(pacbioIndexEnd[i])
		pacbioWrite.WriteString(start + "\t" + end + "\t" + uniquePacbio[i])
	}
}

func genomeFunc(cmd *cobra.Command, args []string) {
	type genomeID struct {
		id string
	}

	type genomeSeq struct {
		seq string
	}

	genomeIDConstruct := []genomeID{}
	genomeSeqConstruct := []genomeSeq{}

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

	genomeIsolate := []string{}

	for i := range genomeSeqConstruct {
		genomeIsolate = append(genomeIsolate, genomeSeqConstruct[i].seq)
	}

	genomeHashes := []string{}
	for i := 0; i <= len(genomeIsolate); i++ {
		for j := 0; j <= len(genomeIsolate[i])-kmerArgs; j++ {
			genomeHashes = append(genomeHashes, genomeIsolate[i][j:j+kmerArgs])
		}
	}

	uniqueGenome, _ := uniqueHash(genomeIsolate)

	genomeFile, err := os.Create("pacbiohashes.txt")
	if err != nil {
		log.Fatal(err)
		defer genomeFile.Close()
	}
	for i := range uniqueGenome {
		genomeFile.WriteString(uniqueGenome[i] + "\n")
	}

	genomeIndexStart := []int{}
	genomeIndexEnd := []int{}

	for i := range uniqueGenome {
		for j := range genomeIsolate {
			start := strings.Index(genomeIsolate[j], uniqueGenome[i])
			end := start + len(uniqueGenome[i])
			genomeIndexStart = append(genomeIndexStart, start)
			genomeIndexEnd = append(genomeIndexEnd, end)
		}
	}

	genomeWrite, err := os.Create("genomeKmerOrigin.txt")
	if err != nil {
		log.Fatal(err)
	}

	for i := range genomeIndexStart {
		start := strconv.Itoa(genomeIndexStart[i])
		end := strconv.Itoa(genomeIndexEnd[i])
		genomeWrite.WriteString(start + "\t" + end + "\t" + uniqueGenome[i])
	}
}

func illuminaFunc(cmd *cobra.Command, args []string) {
	type illuminaID struct {
		id string
	}

	type illuminaSeq struct {
		seq string
	}

	illuminaIDConstruct := []illuminaID{}
	illuminaSeqConstruct := []illuminaSeq{}

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

	illuminaIsolate := []string{}

	for i := range illuminaSeqConstruct {
		illuminaIsolate = append(illuminaIsolate, illuminaSeqConstruct[i].seq)
	}

	illuminaHashes := []string{}

	for i := 0; i <= len(illuminaIsolate); i++ {
		for j := 0; j <= len(illuminaIsolate[i])-kmerArgs; j++ {
			illuminaHashes = append(illuminaHashes, illuminaIsolate[i][j:j+kmerArgs])
		}
	}

	uniqueillumina, _ := uniqueHash(illuminaIsolate)

	illuminaFile, err := os.Create("pacbiohashes.txt")
	if err != nil {
		log.Fatal(err)
		defer illuminaFile.Close()
	}
	for i := range uniqueillumina {
		illuminaFile.WriteString(uniqueillumina[i] + "\n")
	}

	illuminaIndexStart := []int{}
	illuminaIndexEnd := []int{}

	for i := range uniqueillumina {
		for j := range illuminaIsolate {
			start := strings.Index(illuminaIsolate[j], uniqueillumina[i])
			end := start + len(uniqueillumina[i])
			illuminaIndexStart = append(illuminaIndexStart, start)
			illuminaIndexEnd = append(illuminaIndexEnd, end)
		}
	}

	illuminaWrite, err := os.Create("illuminaKmerOrigin.txt")
	if err != nil {
		log.Fatal(err)
	}

	for i := range illuminaIndexStart {
		start := strconv.Itoa(illuminaIndexStart[i])
		end := strconv.Itoa(illuminaIndexEnd[i])
		illuminaWrite.WriteString(start + "\t" + end + "\t" + uniqueillumina[i])
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
