package hw3
import common.execCmd
import java.io.File

import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SAMRecordIterator
import htsjdk.samtools.SamReader
import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.reference.FastaSequenceFile


val minimapPath = "/Johnny/m_vinnichenko/minimap2/minimap2"
val datasetDir = "/Johnny/students/NGS/data/6"
val outDir = "hw3_output"
val readLen = 100

fun Boolean.toInt() = if (this) 1 else 0

class Dataset(val reference: String, val reads: List<String>)

fun runMinimap2(reference: String, reads: List<String>, outFile: String) {
    val cmd = "$minimapPath -ax sr $reference ${reads.joinToString(" ")} > $outFile"
    execCmd(cmd)
}

fun runQuake(reads: List<String>): List<String> {
    val readsCorrected = reads.map { it.replace(".fastq", ".cor.fastq") }
    if (readsCorrected.all { File(it).exists() }) {
        return readsCorrected
    }

    val readNamesFile = "$outDir/reads_file"
    File(readNamesFile).printWriter().use {
        it.println(reads.joinToString(" "))
    }
    val cmd = "python2 $datasetDir/Quake/bin/quake.py -k 15 -f $readNamesFile -p 16"
    execCmd(cmd)

    return readsCorrected
}

fun main(args: Array<String>) {
    println("mainStarted")
    val dataset10K = Dataset(
            reference = "$datasetDir/MG1655-K12.first10K.fasta",
            reads = listOf("$datasetDir/ecoli_10K_err_1.fastq", "$datasetDir/ecoli_10K_err_2.fastq")
    )

    val dataset400K = Dataset(
            reference = "$datasetDir/MG1655-K12.first400K.fasta",
            reads = listOf("$datasetDir/ecoli_400K_err_1.fastq", "$datasetDir/ecoli_400K_err_2.fastq")
    )

    val dataset = dataset400K

    val reads: List<String> = dataset.reads.map { readsFile ->
        val newName = "$outDir/${File(readsFile).getName()}"
        if (!File(newName).exists()) {
            File(readsFile).copyTo(File(newName))
        }
        newName
    }

    val readsCorrected = runQuake(reads)
    val samFile = "$outDir/aln.sam"
    val samCorrectedFile = "$outDir/aln_corrected.sam"

    val referenceSeq: String = FastaSequenceFile(File(dataset.reference), true).use {
        it.nextSequence().baseString
    }

    runMinimap2(dataset.reference, reads, samFile)
    runMinimap2(dataset.reference, readsCorrected, samCorrectedFile)

    val samReaders: Array<SamReader> = arrayOf(
            SamReaderFactory.makeDefault().open(File(samFile)),
            SamReaderFactory.makeDefault().open(File(samCorrectedFile))
    )

    val samIterators: Array<SAMRecordIterator> = samReaders.map { it.iterator() }.toTypedArray()
    val samMaps: Array<MutableMap<String, SAMRecord>> = arrayOf(
            mutableMapOf(),
            mutableMapOf()
    )

    val errorsCount = Array<IntArray>(2, { IntArray(2) })

    while (samIterators.any { it.hasNext() }) {
        for (i in samReaders.indices) {
            if (samIterators[i].hasNext()) {
                val samRecord: SAMRecord = samIterators[i].next()
                if ((samRecord.alignmentBlocks.firstOrNull()?.getLength() ?: 0) < 0.5 * readLen) {
                    continue
                }
                val name: String = samRecord.readName + "/" + samRecord.firstOfPairFlag
                samMaps[i].put(name, samRecord)
                if (samMaps.all { it.containsKey(name) }) {
                    val samRecords: Array<SAMRecord> = arrayOf(samMaps[0][name]!!, samMaps[1][name]!!)
                    if (samRecords[0].alignmentStart == samRecords[1].alignmentStart) {
                        val readSeqs: Array<String> = samRecords.map { it.readString }.toTypedArray<String>()
                        val firstBlocks = samRecords.map { it.alignmentBlocks.first() }
                        val len = samRecords.map {
                            it.alignmentBlocks.map { block ->
                                block.getLength()
                            }.sum()
                        }.min()!!
                        val readStarts = firstBlocks.map { it.getReadStart() - 1 }
                        require(firstBlocks[0].getReferenceStart() == samRecords[0].alignmentStart )
                        require(firstBlocks[1].getReferenceStart() == samRecords[0].alignmentStart )

                        var errors0 = 0
                        var errors1 = 0

                        for (p in 0 until len) {
                            val pos = p + samRecords[0].alignmentStart - 1
                            if (pos >= referenceSeq.length) {
                                break
                            }
                            val error0 = (readSeqs[0][p + readStarts[0]] != referenceSeq[pos]).toInt()
                            val error1 = (readSeqs[1][p + readStarts[1]] != referenceSeq[pos]).toInt()
                            errors0 += error0
                            errors1 += error1
                            errorsCount[error0][error1]++
                        }
                        val threshold = 0.2 * readLen
                        val maxErrors = kotlin.math.max(errors0, errors1)
                        require(maxErrors < threshold)
                    }
                    samMaps.forEach {
                        it.remove(name)
                    }
                }
            }
        }
    }

    println()
    println("false positives ${errorsCount[1][1]}")
    println("false negatives ${errorsCount[0][1]}")
    println("true  negatives ${errorsCount[1][0]}")
}
