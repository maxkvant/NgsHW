package hw4

import common.execCmd
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReader
import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.fastq.FastqReader
import htsjdk.samtools.reference.FastaSequenceFile
import java.io.File
import java.util.*
import kotlin.math.max
import kotlin.math.min

const val outDir = "hw4_output"
const val insertSizeThreshold = 2000

class Dataset(val reference: String, val leftReads: String, val rightReads: String)

fun runPairedBowtie2(dataset: Dataset, samFileName: String) {
    val indexFileName = "index"
    val exitCode1 = execCmd("bowtie2-build ${dataset.reference} $indexFileName ")
    require(exitCode1 == 0)
    val cmd =
            "bowtie2 -k 20 -x $indexFileName --very-sensitive --threads 16 " +
            "-1 ${dataset.leftReads} -2 ${dataset.rightReads} " +
            "> $samFileName"

    val exitCode2 = execCmd(cmd)
    //require(exitCode2 == 0)
}

fun insertSize(left: SAMRecord, right: SAMRecord): Int {
    return max(left.alignmentEnd, right.alignmentEnd) - min(left.alignmentStart, right.alignmentStart)
}

fun forSamRecords(samFile: String, f: (SAMRecord, SAMRecord) -> Unit) {
    val samRecordIterator = SamReaderFactory.makeDefault().open(File(samFile)).iterator()
    while (samRecordIterator.hasNext()) {
        val left = samRecordIterator.next()
        if (!samRecordIterator.hasNext()) {
            continue
        }
        val right = samRecordIterator.next()
        if (insertSize(left, right) < insertSizeThreshold) {
            f(left, right)
        }
    }
}

fun fixName(name: String): String {
    return if (name.endsWith("/1") or name.endsWith("/2")) {
        name.substring(0, name.length - 2)
    } else {
        name
    }
}


fun runTask1(dataset: Dataset, datasetName: String) {
    println("task1")

    val samFile = "$outDir/${datasetName}_aln.sam"
    runPairedBowtie2(dataset, samFile)

    val insertSizes: MutableMap<Int, Int> = TreeMap()

    val referenceStr: String = FastaSequenceFile(File(dataset.reference), false).use {
        it.nextSequence().baseString
    }

    val fastqNames = mutableSetOf<String>()

    for (fastqFile in listOf(dataset.leftReads, dataset.rightReads)) {
        FastqReader(File(fastqFile)).use {
            it.forEach { fastqRecord ->
                fastqNames.add(fixName(fastqRecord.readName))
            }
        }
    }

    val n = referenceStr.length
    val coverages = Array(n, { 0 })

    forSamRecords(samFile, { left, right ->
        val insertSize = insertSize(left, right)
        left.readName = fixName(left.readName)
        right.readName = fixName(right.readName)
        require(left.readName == right.readName)
        insertSizes[insertSize] = (insertSizes[insertSize] ?: 0) + 1
    })

    insertSizes.forEach {insertSize, count ->
        println("insertSize=$insertSize count=$count")
    }
}

fun main(args: Array<String>) {
    val dataDir = "/Johnny/students/NGS/data/4/E.coli"
    val testDataset = Dataset(
            "$dataDir/MG1655-K12.first10K.fasta",
            "$dataDir/test_1.fastq",
            "$dataDir/test_2.fastq"
    )

    val runDataset1 = Dataset(
            "$dataDir/MG1655-K12.fasta",
            "$dataDir/ecoli_mda_lane1_left.fastq.00.cor.fastq",
            "$dataDir/ecoli_mda_lane2_left.fastq.00.cor.fastq"
    )

    val runDataset2 = Dataset(
            "$dataDir/MG1655-K12.fasta",
            "$dataDir/s_6_1.fastq",
            "$dataDir/s_6_2.fastq"
    )

    File(outDir).mkdirs()

    runTask1(testDataset, "test")
}