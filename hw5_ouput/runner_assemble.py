#!/usr/bin/python3.5

import subprocess

import os


class Dataset:
    def __init__(self, name, reads, long_reads, reference):
        self.name = str(name)
        self.reads_1, self.reads_2 = reads
        self.long_reads = long_reads
        self.reference = reference
    def __str__(self):
        return str((self.reads_1, self.reads_2, self.long_reads, self.reference))

reads = ("/Johnny/students/NGS/data/7/illumina.100x.1.fq", "/Johnny/students/NGS/data/7/illumina.100x.2.fq")
ref = "/Johnny/students/NGS/data/7/reference.fasta"
datasets = {"10x": Dataset("10x", reads, "/Johnny/students/NGS/data/7/pacbio_10x.fq", ref),
            "20x": Dataset("20x", reads, "/Johnny/students/NGS/data/7/pacbio_20x.fq", ref),
            "40x": Dataset("40x", reads, "/Johnny/students/NGS/data/7/pacbio_40x.fq", ref)}

class Assembler:
    def __init__(self, dataset, outdir):
        self.dataset = dataset
        self.big_outdir = outdir
        self.outdir = outdir + "/run"

    def getArgs(self):
        return []

    def run(self):
        args = self.getArgs()
        print(args)
        subprocess.call(args)
        print(self.getContigs())

    def getContigs(self):
        pass


class Canu(Assembler):
    def __init__(self, dataset, outdir):
        super().__init__(dataset, outdir + "/canu")

    def getArgs(self):
        canu_path = "/home/m_vinnichenko/intership/canu-1.5/Linux-amd64/bin/canu"
        args = [canu_path]
        args += ["-p", self.dataset.name]
        args += ["maxThreads=24"]
        args += ["-d", self.outdir]
        args += ["useGrid=false"]
        args += ["genomeSize=5m"]
        args += ["-nanopore-raw", self.dataset.long_reads]
        return args

    def getContigs(self):
        return "{}/{}.contigs.fasta".format(self.outdir, self.dataset.name)

class Spades(Assembler):
    def __init__(self, dataset, outdir, pref):
        super().__init__(dataset, outdir + "/" + pref + "spades")

    def getArgs(self):
        spades_path = "/home/m_vinnichenko/intership/algorithmic-biology/assembler/spades.py"
        args = [spades_path]
        args += ["-o", self.outdir]
        args += ["-1", self.dataset.reads_1]
        args += ["-2", self.dataset.reads_2]
        args += ['-t', '12']
        return args

    def getContigs(self):
        return self.outdir + "/contigs.fasta"

class HybridSpades(Spades):
    def __init__(self, dataset, outdir):
        super().__init__(dataset, outdir, "hybrid_")

    def getArgs(self):
        return super().getArgs() + ['--nanopore', self.dataset.long_reads]


def run():
        try:
            spades = Spades(datasets["10x"], "10x", "")
            canu = Canu(datasets["40x"], "40x")
            hybridSpades = HybridSpades(datasets["10x"], "10x")
            hybridSpades20x = HybridSpades(datasets["20x"], "20x")
            assemblers = [spades, hybridSpades, hybridSpades20x, canu]

            for assembler in assemblers:
                assembler.run()

        except Exception as e:
            print("! ERROR " + dataset.name)
            print(e)

if __name__ == '__main__':
    run()
