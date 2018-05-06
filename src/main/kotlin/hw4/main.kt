package hw4

import common.execCmd
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import java.io.File
import kotlin.math.max
import kotlin.math.min

const val insertSizeThreshold = 2000

fun runBowtie2(dataset: Dataset, samFile: String) {
    require(samFile.endsWith(".sam"))
    val samFileName = samFile.substring(0, samFile.length - 4)
    val indexFile = "${samFileName}_index"
    val exitCode1 = execCmd("bowtie2-build ${dataset.reference} $indexFile ")
    require(exitCode1 == 0)

    val cmd = if (dataset.reads.size == 2) {
        "bowtie2 -k 20 -x $indexFile --very-sensitive --threads 16 " +
                "-1 ${dataset.reads[0]} -2 ${dataset.reads[1]} " +
                "> $samFile"
    } else {
        "bowtie2 -k 20 -x $indexFile --very-sensitive --threads 16 " +
                "-U ${dataset.reads[0]} " +
                "> $samFile"

    }
    val exitCode2 = execCmd(cmd)
    println("exitCode = $exitCode2")
}

fun runBwa(dataset: Dataset, samFile: String) {
    require(samFile.endsWith(".sam"))
    val samFileName = samFile.substring(0, samFile.length - 4)
    val exitCode1 = execCmd("bwa index ${dataset.reference}")
    require(exitCode1 == 0)

    val cmd = "bwa mem ${dataset.reference} ${dataset.reads.joinToString(" ")}"
    val exitCode2 = execCmd(cmd)
    println("exitCode = $exitCode2")
}


fun insertSize(left: SAMRecord, right: SAMRecord): Int {
    return max(left.alignmentEnd, right.alignmentEnd) - min(left.alignmentStart, right.alignmentStart)
}

fun fixName(name: String): String {
    return if (name.endsWith("/1") or name.endsWith("/2")) {
        name.substring(0, name.length - 2)
    } else {
        name
    }
}

fun forSamRecordsPaired(samFile: String, f: (SAMRecord, SAMRecord) -> Unit) {
    val samRecordIterator = SamReaderFactory.makeDefault().open(File(samFile)).iterator()
    var iter = 0

    while (samRecordIterator.hasNext()) {
        iter += 1
        if (iter % 1000000 == 0) {
            println("iteration $iter")
        }

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

fun forSamRecordsSingle(samFile: String, paired: Boolean, f: (SAMRecord) -> Unit) {
    if (paired) {
        forSamRecordsPaired(samFile, { left, right ->
            f(left)
            f(right)
        })
    } else {
        val samRecordIterator = SamReaderFactory.makeDefault().open(File(samFile)).iterator()
        var iter = 0
        samRecordIterator.forEach {
            iter += 1
            if (iter % 1000000 == 0) {
                println("iteration $iter")
            }
            f(it)
        }
    }
}

fun runTask1(dataset: Dataset, datasetName: String, directory: String) {
    val samFile = "$directory/${datasetName}_aln.sam"
    if (!File(samFile).exists()) {
        runBowtie2(dataset, samFile)
    }
    val processor = DatasetProcessor(
            dataset = dataset,
            datasetName = datasetName,
            samFile = samFile,
            directory = directory
    )

    processor.plotCoverage()

    val tasks = listOf(
            { processor.countAlignedReads() },
            { processor.plotInsertSizes() },
            { processor.countSubstitutionTable() }
    )
    tasks.parallelStream().forEach {
        it.invoke()
    }
}

fun main(args: Array<String>) {
    val outDir = "hw4_output"

    val dataDir = "/Johnny/students/NGS/data/4/E.coli"
    val testDataset = Dataset(
            "$dataDir/MG1655-K12.first10K.fasta",
            listOf("$dataDir/test_1.fastq",
                    "$dataDir/test_2.fastq")
    )

    val runDataset1 = Dataset(
            "$dataDir/MG1655-K12.fasta",
            listOf("$dataDir/ecoli_mda_lane1_left.fastq.00.cor.fastq",
                    "$dataDir/ecoli_mda_lane1_right.fastq.00.cor.fastq")
    )

    val runDataset2 = Dataset(
            "$dataDir/MG1655-K12.fasta",
            listOf("$dataDir/s_6_1.fastq",
                    "$dataDir/s_6_2.fastq")
    )

    File(outDir).mkdirs()

    val tasks = listOf(
        { runTask1(testDataset, "test", outDir) },
        { runTask1(runDataset1, "e_coli_mda_lane", outDir) },
        { runTask1(runDataset2, "e_coli_s_6", outDir) }
    )
    tasks[0].invoke()
    tasks.parallelStream().forEach {
        it.invoke()
    }
}