package hw4

import common.execCmd
import htsjdk.samtools.SAMRecord
import kotlin.math.max
import kotlin.math.min

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
    val exitCode1 = execCmd("bwa index ${dataset.reference}")
    require(exitCode1 == 0)

    val cmd = "bwa mem -t 16 ${dataset.reference} ${dataset.reads.joinToString(" ")} > $samFile"
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

enum class Aligner {
    BWA, Bowtie2
}
