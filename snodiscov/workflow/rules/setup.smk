# =============================================================================
#  Setup: everything needing internet access. All rules are local rules, so
#  they run in the main Snakemake process (login / interactive node), never on
#  compute nodes. Run them ahead of time with `snodiscov setup`.
# =============================================================================
_REMOTE = {k: s for k, s in config["resources"].items() if "path" not in s and "package" not in s}
_FNAME2KEY = {os.path.basename(RES[k]): k for k in _REMOTE}


localrules:
    fetch_resource,
    install_coco,
    install_pairedBamToBed12,
    setup,


rule fetch_resource:
    """Download a reference file from Zenodo (MD5-verified)."""
    output:
        os.path.join(RESDIR, "data", "{fname}"),
    wildcard_constraints:
        fname=_alt(_FNAME2KEY),
    params:
        spec=lambda wc: config["resources"][_FNAME2KEY[wc.fname]],
    retries: 2
    run:
        from snodiscov.fetch import fetch_resource

        key = _FNAME2KEY[wildcards.fname]
        try:
            fetch_resource(params.spec, output[0])
        except Exception as e:
            raise RuntimeError(
                "resource '" + key + "' " + str(params.spec) + ": " + str(e)
                + "\n(check resources: " + key + " in the run config.yaml)"
            ) from e


rule install_coco:
    """CoCo (Deschamps-Francoeur et al. 2019) pinned at a fixed git ref."""
    output:
        touch(COCO_DONE),
    params:
        repo=config["tools"]["coco"]["repo"],
        ref=config["tools"]["coco"]["ref"],
        dest=os.path.join(TOOLS, "coco"),
    log:
        os.path.join(TOOLS, "logs", "install_coco.log"),
    conda:
        ENVS + "/git.yaml"
    shell:
        """
        (rm -rf {params.dest}
         git clone {params.repo} {params.dest}
         git -C {params.dest} checkout {params.ref}
         git -C {params.dest} rev-parse HEAD > {params.dest}/COMMIT) &> {log}
        """


rule install_pairedBamToBed12:
    output:
        touch(PB_DONE),
    params:
        repo=config["tools"]["pairedBamToBed12"]["repo"],
        ref=config["tools"]["pairedBamToBed12"]["ref"],
        dest=os.path.join(TOOLS, "pairedBamToBed12"),
    log:
        os.path.join(TOOLS, "logs", "install_pairedBamToBed12.log"),
    conda:
        ENVS + "/git.yaml"
    shell:
        """
        (rm -rf {params.dest}
         git clone {params.repo} {params.dest}
         git -C {params.dest} checkout {params.ref}
         git -C {params.dest} rev-parse HEAD > {params.dest}/COMMIT
         make -C {params.dest}
         test -x {params.dest}/bin/pairedBamToBed12) &> {log}
        """


rule setup:
    """Target used by `snodiscov setup`: all downloads + tool installs."""
    input:
        [RES[k] for k in _REMOTE],
        COCO_DONE,
        PB_DONE if COVERAGE else [],
