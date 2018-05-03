package hw4

import common.execCmd
import com.github.sh0nk.matplotlib4j.Plot
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.fastq.FastqReader
import htsjdk.samtools.reference.FastaSequenceFile
import java.io.File
import java.io.PrintWriter
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
    print("exitCode = $exitCode2")
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

const val nucleotides = "acgt"

fun nucleotideCode(c: Char): Int = "${nucleotides}n".indexOf(c.toLowerCase())

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
    val coverage = Array(n, { 0 })
    val coverageOtherVariants = Array(n, { 0 })
    val substitutionCount: Array<Array<Int>> = Array(4, { Array(4, { 0 }) })

    val allReads = fastqNames.size

    forSamRecords(samFile, { left, right ->
        val insertSize = insertSize(left, right)
        left.readName = fixName(left.readName)
        right.readName = fixName(right.readName)
        require(left.readName == right.readName)
        insertSizes[insertSize] = (insertSizes[insertSize] ?: 0) + 1
        fastqNames.remove(left.readName)

        for (samRead in arrayOf(left, right)) {
            for (alignmentBlock in samRead.alignmentBlocks) {
                for (pos in 0 until alignmentBlock.length) {
                    val posRead = alignmentBlock.readStart - 1
                    val posReference = alignmentBlock.referenceStart - 1
                    val cRead = nucleotideCode(samRead.readString[posRead])
                    val cRef = nucleotideCode(referenceStr[posReference])
                    coverage[posReference] += 1
                    if (cRead != cRead) {
                        coverageOtherVariants[posReference] += 1
                    }
                    if (cRead >= 4 || cRef >= 4) {
                        continue
                    }
                    substitutionCount[cRef][cRead] += 1
                }
            }
        }
    })

    forSamRecords(samFile, { left, right ->
        val insertSize = insertSize(left, right)
        left.readName = fixName(left.readName)
        right.readName = fixName(right.readName)
        require(left.readName == right.readName)
        insertSizes[insertSize] = (insertSizes[insertSize] ?: 0) + 1
        fastqNames.remove(left.readName)

        for (samRead in arrayOf(left, right)) {
            for (alignmentBlock in samRead.alignmentBlocks) {
                for (pos in 0 until alignmentBlock.length) {
                    val posRead = alignmentBlock.readStart - 1
                    val posReference = alignmentBlock.referenceStart - 1
                    val cRead = nucleotideCode(samRead.readString[posRead])
                    val cRef = nucleotideCode(referenceStr[posReference])
                    coverage[posReference] += 1
                    if (cRef != cRead) {
                        coverageOtherVariants[posReference] += 1
                    }
                }
            }
        }
    })

    forSamRecords(samFile, { left, right ->
        for (samRead in arrayOf(left, right)) {
            for (alignmentBlock in samRead.alignmentBlocks) {
                for (pos in 0 until alignmentBlock.length) {
                    val posRead = alignmentBlock.readStart - 1
                    val posReference = alignmentBlock.referenceStart - 1
                    val cRead = nucleotideCode(samRead.readString[posRead])
                    val cRef = nucleotideCode(referenceStr[posReference])
                    val covOther = coverageOtherVariants[posReference]
                    val isSnp = covOther >= 3 && covOther > coverage[posReference] * 0.3
                    if (isSnp || cRead >= 4 || cRef >= 4 || cRead == cRef) {
                        continue
                    }
                    substitutionCount[cRef][cRead] += 1
                }
            }
        }
    })

    val lostReads = fastqNames.size
    val alignedReads = allReads - lostReads

    val outFile = "$outDir/${datasetName}_out.txt"
    val printer = PrintWriter(File(outFile))
    printer.use {
        it.println("readsAligned: $alignedReads (${alignedReads / allReads.toDouble()})")
        it.println(nucleotides.toCharArray().joinToString(" ", "  "))
        for (i in substitutionCount.indices) {
            substitutionCount[i][i] = 0
            it.print("${nucleotides[i]} ${substitutionCount[i].joinToString(" ")}")
            it.println()
        }
    }

    fun plotInsertSize() {
        val plt: Plot = Plot.create()
        plt.ylabel("# reads")
        plt.xlabel("insertSize")
        plt.plot().add(insertSizes.map { it.key },
                insertSizes.map { it.value },
                "r+")

        plt.savefig("$outDir/${datasetName}_insert_size.png").dpi(200.0)
        plt.executeSilently()
    }

    fun plotCoverage() {
        val coverageAverages = mutableListOf<Double>()
        val positions = mutableListOf<Double>()
        val stepLen = 50
        val windowSize = 1000
        for (i in 0 until coverage.size - windowSize step stepLen) {
            val average: Double = coverage.sliceArray(i until min(coverage.size, i + windowSize)).average()
            coverageAverages.add(average)
            positions.add(i + windowSize / 2.0)
        }

        val plt: Plot = Plot.create()
        plt.ylabel("coverage")
        plt.xlabel("pos")
        plt.plot().add(positions,
                coverageAverages,
                "r+")

        plt.savefig("$outDir/${datasetName}_coverage.png").dpi(200.0)
        plt.executeSilently()
    }

    plotInsertSize()
    plotCoverage()
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
            "$dataDir/ecoli_mda_lane1_right.fastq.00.cor.fastq"
    )

    val runDataset2 = Dataset(
            "$dataDir/MG1655-K12.fasta",
            "$dataDir/s_6_1.fastq",
            "$dataDir/s_6_2.fastq"
    )

    File(outDir).mkdirs()

    //runTask1(testDataset, "test")

    runTask1(runDataset1, "e_coli_mda_lane")
    runTask1(runDataset2, "e_coli_s_6")
}