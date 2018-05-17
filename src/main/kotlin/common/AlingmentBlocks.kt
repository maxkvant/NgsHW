package common

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

