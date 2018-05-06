package hw4

data class Dataset(val reference: String, val reads: List<String>) {
    val paired
            get(): Boolean = reads.size == 2

    init {
        require(reads.size == 1 || reads.size == 2)
    }
}