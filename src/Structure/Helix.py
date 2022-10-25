from attrs import field, define
import mpmath as mpm

@define
class HelixBase:
    compounds: dict = define(kw_only=True, repr=False)
    calculator: str = define(kw_only=True, repr=False)
