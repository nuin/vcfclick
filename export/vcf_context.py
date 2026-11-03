"""Typed request and collected metadata shared by the VCF export stages."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

Region = tuple[str, int | None, int | None]
InfoColumn = tuple[str, str, int]
SiteCalls = dict[str, list[Any]]
FORMAT_PLACEHOLDER = "\x00FMT\x00"
FORMAT_FIELDS = ("GT", "GQ", "DP", "AD", "FT")


class ExportError(ValueError):
    """A user-facing problem with the export request."""


@dataclass(frozen=True)
class ExportRequest:
    db_name: str
    ingest_id: str | None = None
    samples: list[str] | None = None
    regions: list[Region] | None = None
    where: str | None = None
    pass_only: bool = False
    sites_only: bool = False
    absent_as: str = "ref"
    version: str = ""

    def filters(self) -> dict[str, Any]:
        return {
            "regions": self.regions,
            "where": self.where,
            "pass_only": self.pass_only,
            "samples": self.samples,
            "sites_only": self.sites_only,
        }


@dataclass
class ExportContext:
    request: ExportRequest
    ingest_id: str
    samples: list[str]
    keep_ref: bool
    absent_gt: str
    variant_where: str
    info_cols: list[InfoColumn]
    has_extra: bool
    chrom_stats: list[list[Any]]
    sample_index: dict[str, int] = field(init=False)
    used_info: dict[str, str] = field(default_factory=dict)
    used_extra: set[str] = field(default_factory=set)
    used_filters: set[str] = field(default_factory=set)
    fmt_seen: set[str] = field(default_factory=set)

    def __post_init__(self) -> None:
        self.sample_index = {sample: i for i, sample in enumerate(self.samples)}

    @property
    def select_columns(self) -> list[str]:
        return (
            ["chrom", "pos", "vcf_id", "ref", "alt", "qual", "filter"]
            + [column for column, _, _ in self.info_cols]
            + (["info_extra"] if self.has_extra else [])
        )
