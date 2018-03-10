package hw2.genome_assembler

import htsjdk.samtools.fastq.FastqReader
import htsjdk.samtools.fastq.FastqRecord
import java.io.File

typealias Kmer = String
val K = 55
fun Kmer.reverseComplement(): Kmer {
    val chars = this.toCharArray()
    chars.reverse()
    for (i in chars.indices) {
        val c = chars[i]
        chars[i] = when (c) {
            'a' -> 't'
            't' -> 'a'
            'g' -> 'c'
            'c' -> 'g'
            else -> 'n'
        }
    }
    return String(chars)
}

fun main(args: Array<String>) {
    val startTime = System.nanoTime();

    val pathTest = "/Johnny/data/input/Bacteria/E.coli/K12/is220/cropped/s_6.first10000_1.fastq.gz"

    val file = File(pathTest)

    val kmers: MutableMap<Kmer, Int> = mutableMapOf()

    FastqReader(file).use {
        var reads = 0
        it.iterator().forEach { fastqRecord: FastqRecord ->
            reads++
            val infoCount = 100000
            if (reads % infoCount == 0) {
                println("processed $reads reads")
            }

            val readChars = fastqRecord.readString.toLowerCase().toCharArray()
            for (i in 0 .. readChars.size - K) {
                val kmer1: Kmer = String(readChars, i, K)
                val kmer2 = kmer1.reverseComplement()
                kmers[kmer1] = (kmers[kmer1] ?: 0) + 1
                kmers[kmer2] = (kmers[kmer2] ?: 0) + 1
            }
        }
    }

    print("${System.nanoTime() - startTime}")

}