package common

import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import hw4.insertSize
import java.io.File

const val insertSizeThreshold = 2000

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

