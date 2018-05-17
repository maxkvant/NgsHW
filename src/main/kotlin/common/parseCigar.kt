package common

import htsjdk.samtools.CigarOperator
import htsjdk.samtools.SAMRecord

fun SAMRecord.parseCigar(referenceChars: CharArray): List<AlignmentBlock> {
    var readPos = 0
    var referencePos = this.alignmentStart - 1
    val ans = mutableListOf<AlignmentBlock>()
    val readChars = readString.toCharArray()

    for (cigarElement in this.cigar.cigarElements) {
        val op = cigarElement.operator!!
        when (op) {
            CigarOperator.M, CigarOperator.EQ, CigarOperator.X -> {
                val seq = String(readChars, readPos, cigarElement.length)
                ans.add(Match(seq, referencePos, readPos))
            }

            CigarOperator.I -> {
                val seq = String(readChars, readPos, cigarElement.length)
                ans.add(Insertion(seq, referencePos, readPos))
            }

            CigarOperator.D -> {
                val seq = String(referenceChars, referencePos, cigarElement.length)
                ans.add(Deletion(seq, referencePos, readPos))
            }
            CigarOperator.H, CigarOperator.N, CigarOperator.P, CigarOperator.S -> {}
        }
        if (op.consumesReadBases()) {
            readPos += cigarElement.length
        }
        if (op.consumesReferenceBases()) {
            referencePos += cigarElement.length
        }
    }
    return ans
}