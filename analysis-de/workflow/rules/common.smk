def references(self):
        self.refGz = os.path.join(
            ".".join([
                self.species.capitalize(),
                self.assembly,
                "dna.primary_assembly.fa.gz"
            ])
        )
        self.refFa = self.refGz.rstrip(".gz")
        self.refGz_path = os.path.join(self.refs_dir, self.refGz)
        self.refFa_path = os.path.join(self.refs_dir, self.refFa)
        self.ref_url = os.path.join(
            "rsync://ftp.ensembl.org/ensembl/pub",
            "release-" + str(self.ensembl_release),
            "fasta",
            self.species,
            "dna",
            self.refGz
        )
        self.gtf = os.path.join(
            ".".join([
                self.species.capitalize(),
                self.assembly,
                str(self.ensembl_release),
                "chr.gtf.gz"
            ])
        )
        self.gtf_path = os.path.join(self.refs_dir, self.gtf)
        self.gtf_url = os.path.join(
            "rsync://ftp.ensembl.org/ensembl/pub",
            "release-" + str(self.ensembl_release),
            "gtf",
            self.species,
            self.gtf
        )