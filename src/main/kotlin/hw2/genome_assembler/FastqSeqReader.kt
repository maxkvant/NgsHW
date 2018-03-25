package hw2.genome_assembler

import htsjdk.samtools.fastq.FastqReader
import java.io.File

class FastqSeqReader(private val file: File): SeqReader {
    override fun forStrings(f: (String) -> Unit) {
        FastqReader(file).use {
            it.forEach { fastqRecord ->
                f(fastqRecord.readString)
            }
        }
    }
}