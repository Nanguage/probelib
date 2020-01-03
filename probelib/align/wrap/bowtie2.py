import typing as t
import subprocess as subp


def align_se_sen(fq_path: str,
                 index: str,
                 sam_path: str,
                 threads: int = 10,
                 log: t.Optional[str] = 'bowtie2.log',
                 header: bool = True,
                 ) -> str:
    cmd = ["bowtie2", "-x", index, "-U", fq_path]
    if not header:
        cmd.append("--no-hd")
    cmd += ["-t", "-k", "100", "--very-sensitive-local"]
    cmd += ["-p", str(threads)]
    cmd += ["-S", sam_path]
    cmd = " ".join(cmd)
    if log:
        cmd += f" 2> {log}"
    subp.check_call(cmd, shell=True)
    return sam_path
