package hw2.genome_assembler

import reverseComplement

interface SeqReader: Cloneable {
    fun forStrings(f: (String) -> Unit)

    fun forCharArraysComplement(f: (CharArray) -> Unit) {
        forStrings { str ->
            val seq1 = str.toLowerCase()
            val seq2 = reverseComplement(seq1)

            f(seq1.toCharArray())
            f(seq2.toCharArray())
        }
    }
}