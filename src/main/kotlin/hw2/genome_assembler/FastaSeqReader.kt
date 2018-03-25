package hw2.genome_assembler

import htsjdk.samtools.reference.FastaSequenceFile
import htsjdk.samtools.reference.ReferenceSequence
import java.io.File

class FastaSeqReader(private val file: File): SeqReader {
    override fun forStrings(f: (String) -> Unit) {
        FastaSequenceFile(file, false).use {
            while (true) {
                val sequence = it.nextSequence()
                if (sequence == null) {
                    break
                } else {
                    f(sequence.baseString)
                }
            }
        }
    }
}