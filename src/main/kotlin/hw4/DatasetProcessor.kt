package hw4

import com.github.sh0nk.matplotlib4j.Plot
import htsjdk.samtools.fastq.FastqReader
import htsjdk.samtools.reference.FastaSequenceFile
import java.io.File
import java.io.PrintWriter
import java.util.*
import kotlin.math.min

class DatasetProcessor(
        val dataset: Dataset,
        val datasetName: String,
        val samFile: String,
        val directory: String
) {
    companion object {
        private const val nucleotides = "acgt_"
        fun nucleotideCode(c: Char): Int = "${nucleotides}n".indexOf(c.toLowerCase())
    }

    @Volatile
    private var coverageStats: CoverageStats? = null

    private val referenceStr = FastaSequenceFile(File(dataset.reference), false).use {
        it.nextSequence().baseString
    }

    private val referenceChars: CharArray = referenceStr.toCharArray()
    private val n = referenceStr.length
    private val paired = dataset.reads.size == 2

    private inner class CoverageStats {
        val coverage = Array(n, { 0 })
        val coverageOtherVariants = Array(n, { 0 })

    }

    private fun getCoverageStats(): CoverageStats = synchronized(this) {
        val curCoverageStats = coverageStats
        if (curCoverageStats != null) {
            return curCoverageStats
        } else {
            println("collecting coverage statics")
            val stats = CoverageStats()
            forSamRecordsSingle(samFile, paired, { samRead ->
                for (alignmentBlock in samRead.alignmentBlocks) {
                    for (pos in 0 until alignmentBlock.length) {
                        val posRead = alignmentBlock.readStart - 1
                        val posReference = alignmentBlock.referenceStart - 1
                        val cRead = nucleotideCode(samRead.readString[posRead])
                        val cRef = nucleotideCode(referenceStr[posReference])
                        stats.coverage[posReference] += 1
                        if (cRef != cRead) {
                            stats.coverageOtherVariants[posReference] += 1
                        }
                    }
                }
            })
            coverageStats = stats
            return stats
        }
    }

    fun countSubstitutionTable() {
        val coverageStats = getCoverageStats()
        val substitutionCount: Array<Array<Int>> = Array(nucleotides.length, { Array(nucleotides.length, { 0 }) })
        println("calculating substitution table")
        val indelLen = mutableMapOf<Int, Int>()

        forSamRecordsSingle(samFile, paired, { samRead ->
            var alignedChars = 0
            var commonChars = 0
            for (block in samRead.parseCigar(referenceChars)) {
                when (block) {
                    is Match ->
                        for (pos in block.seq.indices) {
                            val posRead = block.readBegin + pos
                            val posReference = block.referenceBegin + pos
                            require(samRead.readString[posRead] == block.seq[pos])
                            val cRead = nucleotideCode(samRead.readString[posRead])
                            val cRef = nucleotideCode(referenceStr[posReference])
                            coverageStats.coverage[posReference] += 1
                            if (cRead != cRead) {
                                coverageStats.coverageOtherVariants[posReference] += 1
                            }
                            if (cRead >= nucleotides.length || cRef >= nucleotides.length) {
                                continue
                            }
                            substitutionCount[cRef][cRead] += 1
                            alignedChars++
                            if (cRef == cRead) {
                                commonChars++
                            }
                        }
                    is Insertion ->
                        for (pos in block.seq.indices) {
                            val posRead = block.readBegin + pos
                            require(samRead.readString[posRead] == block.seq[pos])
                            val cRead = nucleotideCode(samRead.readString[posRead])
                            val cRef = nucleotideCode('_')
                            if (cRead >= nucleotides.length || cRef >= nucleotides.length) {
                                continue
                            }
                            substitutionCount[cRef][cRead] += 1
                            alignedChars++
                            indelLen[block.seq.length] = (indelLen[block.seq.length] ?: 0) + 1
                        }
                    is Deletion ->
                        for (pos in block.seq.indices) {
                            val posReference = block.referenceBegin + pos
                            require(referenceStr[posReference] == block.seq[pos])
                            val cRead = nucleotideCode('_')
                            val cRef= nucleotideCode(referenceStr[posReference])
                            if (cRead >= nucleotides.length || cRef >= nucleotides.length) {
                                continue
                            }
                            substitutionCount[cRef][cRead] += 1
                            alignedChars++
                            indelLen[-block.seq.length] = (indelLen[-block.seq.length] ?: 0) + 1
                        }
                }
                //require(commonChars > alignedChars * 0.7)
            }
        })

        val outFile = "$directory/${datasetName}_substitutions.txt"
        PrintWriter(File(outFile)).use {
            it.println(nucleotides.toCharArray().joinToString(" ", "  "))
            for (i in substitutionCount.indices) {
                substitutionCount[i][i] = 0
                it.print("${nucleotides[i]} ${substitutionCount[i].joinToString(" ")}")
                it.println()
            }
        }

        val plt: Plot = Plot.create()
        plt.ylabel("count")
        plt.xlabel("indel len")
        plt.plot().add(indelLen.map { it.key }, indelLen.map { it.value }, "r+")
        plt.savefig("$directory/${datasetName}_indel_len.png").dpi(200.0)
        plt.executeSilently()
    }

    fun plotInsertSizes() {
        require(paired)
        println("plotting insert sizes")
        val insertSizes: MutableMap<Int, Int> = TreeMap()
        forSamRecordsPaired(samFile, { left, right ->
            if (left.readName == right.readName) {
                val insertSize = insertSize(left, right)
                insertSizes[insertSize] = (insertSizes[insertSize] ?: 0) + 1
            }
        })

        val plt: Plot = Plot.create()
        plt.ylabel("# reads")
        plt.xlabel("insertSize")
        plt.plot().add(insertSizes.map { it.key }, insertSizes.map { it.value }, "r+")
        plt.savefig("$directory/${datasetName}_insert_size.png").dpi(200.0)
        plt.executeSilently()

        val reads = insertSizes.map { it.value }.sum()
        var curReads = 0
        var intervalBegin = 0
        var intervalEnd = insertSizeThreshold
        insertSizes.forEach { key, value ->
            curReads += value
            if (curReads <= reads * 0.025) {
                intervalBegin = key
            }
            if (reads - curReads >= reads * 0.025) {
                intervalEnd = key
            }
        }

        val plt2 = Plot.create()
        plt2.ylabel("# reads")
        plt2.xlabel("insertSize")
        plt2.plot().add(insertSizes.map { it.key }, insertSizes.map { it.value }, "r+")
        plt2.title("95% confidence interval")
        plt2.xlim(intervalBegin, intervalEnd)
        plt2.savefig("$directory/${datasetName}_insert_size_2.png").dpi(200.0)
        plt2.executeSilently()
    }

    fun plotCoverage() {
        println("plotting coverage")
        val coverageAverages = mutableListOf<Double>()
        val positions = mutableListOf<Double>()
        val stepLen = 50
        val windowSize = 1000
        val coverage = getCoverageStats().coverage

        for (i in 0 until coverage.size - windowSize step stepLen) {
            val average: Double = coverage.sliceArray(i until min(coverage.size, i + windowSize)).average()
            coverageAverages.add(average)
            positions.add(i + windowSize / 2.0)
        }

        val plt: Plot = Plot.create()
        plt.ylabel("coverage")
        plt.xlabel("pos")
        plt.plot().add(positions, coverageAverages, "r+")
        plt.savefig("$directory/${datasetName}_coverage.png").dpi(200.0)
        plt.executeSilently()


    }

    fun countAlignedReads() {
        val fastqNames = mutableSetOf<String>()

        for (fastqFile in dataset.reads) {
            FastqReader(File(fastqFile)).use {
                it.forEach { fastqRecord ->
                    fastqNames.add(fixName(fastqRecord.readName))
                }
            }
        }

        val allReads = fastqNames.size

        forSamRecordsSingle(samFile, paired, {
            fastqNames.remove(fixName(it.readName))
        })

        val lostReads = fastqNames.size
        val alignedReads = allReads - lostReads

        val outFile = "$directory/${datasetName}_read_count.txt"
        PrintWriter(File(outFile)).use {
            it.println("allReads $allReads")
            it.println("alignedReads $alignedReads: ${alignedReads / allReads.toDouble() * 100}%")
            it.println("lostReads $lostReads: ${lostReads / allReads.toDouble() * 100}%")
        }
    }

    fun plotMonomers() {
        fun extendMonomer(str: String, i: Int): Pair<Int, Int> {
            var l = i
            var r = i
            while (l - 1 >= 0 && str[l - 1] == str[i]) {
                l--
            }
            while (r + 1 < str.length && str[r + 1] == str[r]) {
                r++
            }
            return Pair(l, r + 1)
        }

        val lenThreshold = 20
        val minLen = 3
        val monomerCount = Array<MutableMap<Int,Int> >(lenThreshold + 1, { mutableMapOf() })

        forSamRecordsSingle(samFile, paired, { samRead ->
            for (alignmentBlock in samRead.alignmentBlocks) {
                var pos = 0
                while (pos < alignmentBlock.length) {
                    val posRead = alignmentBlock.readStart - 1
                    val posReference = alignmentBlock.referenceStart - 1
                    val cRead = samRead.readString[posRead]
                    val cRef = referenceStr[posReference]
                    if (cRead == cRef) {
                        val refRange = extendMonomer(referenceStr, posReference)
                        val refLen = refRange.second - refRange.first
                        val readRange = extendMonomer(samRead.readString, posRead)
                        val readLen = readRange.second - readRange.first

                        if (minLen <= refLen && refLen <= lenThreshold) {
                            monomerCount[refLen][readLen - readLen] = (monomerCount[refLen][readLen - readLen] ?: 0) + 1
                        }
                        pos = refRange.second
                    } else {
                        pos += 1
                    }
                }
            }
        })
    }
}