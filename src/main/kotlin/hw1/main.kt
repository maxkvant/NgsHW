package hw1

import com.github.sh0nk.matplotlib4j.Plot
import htsjdk.samtools.fastq.FastqReader
import htsjdk.samtools.fastq.FastqRecord
import java.io.File
import java.lang.Math

fun Int.fromPhredScore(): Double = Math.pow(10.0, -this.toDouble() / 10.0)

fun Double.toPhredScore() = Math.round(-10.0 * Math.log10(this))

fun main(ars: Array<String>) {
    val file = File("/Johnny/data/input/Bacteria/E.coli/K12/is220/cropped/s_6.first10000_1.fastq.gz")
    FastqReader(file).use {
        val readLen = 300
        val qualitySum: Array<Double> = Array(readLen, { 0.0 })
        val readCount: Array<Int> = Array(readLen, {0})
        val gcBuckets = 50
        val gcContents = Array<Int>(gcBuckets + 1, {0})

        var reads = 0
        it.iterator().forEach { fastqRecord: FastqRecord ->
            reads++
            val readString = fastqRecord.readString
            val quality = fastqRecord.baseQualities

            require(readString.length == quality.size)

            var gc = 0
            var nucleotides = 0
            for (i in readString.indices) {
                readCount[i]++
                qualitySum[i] += quality[i].toInt().fromPhredScore()
                if (quality[i] < 4) {
                    continue
                }
                val nucleotide = readString[i].toLowerCase()
                if (nucleotide == 'g' || nucleotide == 'c') {
                    gc++
                }
                if (nucleotide != 'n') {
                    nucleotides++
                }
            }
            if (nucleotides != 0) {
                val gcContent = gc.toDouble() / nucleotides
                val gcBucket = Math.floor(gcContent * gcBuckets).toInt()
                gcContents[gcBucket]++
            }
        }

        val plt: Plot = Plot.create()

        plt.plot().add(gcContents.indices.map { i -> (i + 0.5) / gcBuckets}, gcContents.map(Int::toDouble))
        plt.show()

        for (i in gcContents.indices) {
            println("${i.toDouble() / gcBuckets} - ${(i + 1).toDouble() / gcBuckets} :  ${gcContents[i]}")
        }

        for (i in qualitySum.indices) {
            if (readCount[i] == 0) {
                continue
            }
            println("$i ${(qualitySum[i] / readCount[i]).toPhredScore()}")
        }

    }
}
