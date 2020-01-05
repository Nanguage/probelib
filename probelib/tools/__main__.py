import fire

from .candidates import candidates
from .avoid_otp import avoid_otp

CLI = {
    'candidates': candidates,
    'avoid_otp': avoid_otp,
}

fire.Fire(CLI)
