package hw2.genome_assembler

import htsjdk.samtools.fastq.FastqReader
import java.io.File

class FastqSeqReader(private val file: File): SeqReader {
    override fun forStrings(f: (String) -> Unit) {
        FastqReader(file).use {
            var reads = 0
            it.forEach { fastqRecord ->
                reads++
                val infoCount = 100000
                if (reads % infoCount == 0) {
                    println("processed $reads reads")
                }
                f(fastqRecord.readString)
            }
        }
    }
}