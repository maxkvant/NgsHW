package hw4

import htsjdk.samtools.CigarOperator
import htsjdk.samtools.SAMRecord

interface AlignmentBlock {
    val seq: String
    val referenceBegin: Int
    val readBegin: Int
}

data class Insertion(
        override val seq: String,
        override val referenceBegin: Int,
        override val readBegin: Int
): AlignmentBlock


data class Deletion(
        override val seq: String,
        override val referenceBegin: Int,
        override val readBegin: Int
): AlignmentBlock

data class Match(
        override val seq: String,
        override val referenceBegin: Int,
        override val readBegin: Int
): AlignmentBlock

fun SAMRecord.parseCigar(referenceChars: CharArray): List<AlignmentBlock> {
    var readPos = 0
    var referencePos = this.alignmentStart - 1
    val ans = mutableListOf<AlignmentBlock>()
    val readChars = readString.toCharArray()

    for (cigarElement in this.cigar.cigarElements) {
        val op = cigarElement.operator!!
        when (op) {
            CigarOperator.M, CigarOperator.EQ, CigarOperator.X  -> {
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
            CigarOperator.H, CigarOperator.N, CigarOperator.P, CigarOperator.S-> {}
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