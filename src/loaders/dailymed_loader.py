"""
Drug label loader — Gate 2 Pediatric Safety.
Primary source: openFDA drug/label.json API (free, no registration, returns structured JSON).
Fallback: DailyMed SPL XML endpoint.

The openFDA label API returns the same structured sections as DailyMed
(boxed_warning, contraindications, pediatric_use) in clean JSON —
avoiding the DailyMed XML parsing complexity entirely.
"""

import re
import logging
import requests
from .cache import file_cache

logger = logging.getLogger(__name__)
OPENFDA_LABEL_URL = "https://api.fda.gov/drug/label.json"
DAILYMED_URL = "https://dailymed.nlm.nih.gov/dailymed/services/v2"
TIMEOUT = 30

_PEDIATRIC_CONTRA_PATTERNS = [
    r"contraindicated\s+in\s+(pediatric|paediatric|children|infants?|neonates?|patients?\s+under\s+\d+)",
    r"not\s+(recommended|indicated)\s+for\s+(use\s+in\s+)?(pediatric|paediatric|children|infants?)",
    r"safety\s+and\s+efficacy\s+(have\s+not\s+been\s+established|is\s+not\s+established)\s+in\s+(pediatric|children)",
    r"not\s+approved\s+for\s+use\s+in\s+(pediatric|paediatric|children)",
]


@file_cache("dailymed")
def get_label_text(drug_name: str) -> dict:
    """
    Fetch and parse a drug label for pediatric safety scoring.
    Returns dict with: boxed_warning, contraindications, pediatric_use,
    has_pediatric_bb_warning, has_pediatric_contraindication,
    min_age_years, pediatric_approved.
    """
    # Try openFDA drug label API first (JSON, structured sections)
    label = _fetch_openfda_label(drug_name)
    if label and (label.get("pediatric_use") or label.get("contraindications") or label.get("boxed_warning")):
        return label

    # Fallback: DailyMed XML
    return _fetch_dailymed_xml_label(drug_name)


def _fetch_openfda_label(drug_name: str) -> dict:
    """Fetch label sections from openFDA drug/label.json API."""
    for search_field in [
        f'openfda.generic_name:"{drug_name}"',
        f'openfda.brand_name:"{drug_name}"',
        f'openfda.substance_name:"{drug_name}"',
    ]:
        try:
            resp = requests.get(
                OPENFDA_LABEL_URL,
                params={"search": search_field, "limit": 1},
                timeout=TIMEOUT,
            )
            if resp.status_code == 404:
                continue
            resp.raise_for_status()
            results = resp.json().get("results", [])
            if not results:
                continue
            label = results[0]
            return _parse_openfda_label(label)
        except Exception as e:
            logger.debug("openFDA label search '%s': %s", search_field, e)

    return _empty_label()


def _parse_openfda_label(label: dict) -> dict:
    """Extract pediatric fields from openFDA label result."""
    # openFDA returns section text as lists of strings
    def get_text(key: str) -> str:
        val = label.get(key, [])
        if isinstance(val, list):
            return " ".join(val)
        return str(val) if val else ""

    boxed_warning = get_text("boxed_warning") + " " + get_text("warnings_and_precautions")
    contraindications = get_text("contraindications")
    pediatric_use = get_text("pediatric_use") + " " + get_text("use_in_specific_populations")

    has_bb_warning = _has_pediatric_bb_content(boxed_warning)
    has_contra = _has_pediatric_contraindication(contraindications + " " + pediatric_use)
    min_age = _extract_min_age(pediatric_use)
    approved = _is_pediatric_approved(pediatric_use)

    return {
        "boxed_warning": boxed_warning.strip()[:500],
        "contraindications": contraindications.strip()[:500],
        "pediatric_use": pediatric_use.strip()[:1000],
        "has_pediatric_bb_warning": has_bb_warning,
        "has_pediatric_contraindication": has_contra,
        "min_age_years": min_age,
        "pediatric_approved": approved,
    }


def _fetch_dailymed_xml_label(drug_name: str) -> dict:
    """Fallback: fetch SPL XML from DailyMed and parse pediatric sections."""
    try:
        # Step 1: resolve to setid
        resp = requests.get(
            f"{DAILYMED_URL}/spls.json",
            params={"drug_name": drug_name, "limit": 3},
            timeout=TIMEOUT,
        )
        resp.raise_for_status()
        spls = resp.json().get("data", [])
        if not spls:
            return _empty_label()
        setid = spls[0].get("setid", "")
        if not setid:
            return _empty_label()

        # Step 2: fetch XML
        xml_resp = requests.get(
            f"{DAILYMED_URL}/spls/{setid}.xml",
            timeout=TIMEOUT,
        )
        if xml_resp.status_code != 200:
            return _empty_label()

        return _parse_spl_xml(xml_resp.text)
    except Exception as e:
        logger.debug("DailyMed XML fallback error for %s: %s", drug_name, e)
        return _empty_label()


def _parse_spl_xml(xml_text: str) -> dict:
    """Parse SPL HL7 CDA XML to extract pediatric-relevant sections."""
    import xml.etree.ElementTree as ET
    # SPL section codes
    SECTION_CODES = {
        "34066-1": "boxed_warning",
        "34070-3": "contraindications",
        "34081-0": "pediatric_use",
        "43685-7": "warnings",
    }
    sections = {k: "" for k in SECTION_CODES.values()}

    try:
        # Remove namespace prefixes for simple parsing
        xml_clean = re.sub(r'\sxmlns[^"]*"[^"]*"', '', xml_text)
        xml_clean = re.sub(r'<([a-zA-Z]+):([a-zA-Z]+)', r'<\2', xml_clean)
        xml_clean = re.sub(r'</([a-zA-Z]+):([a-zA-Z]+)', r'</\2', xml_clean)

        root = ET.fromstring(xml_clean)

        for section in root.iter("section"):
            code_elem = section.find("code")
            if code_elem is None:
                continue
            code = code_elem.get("code", "")
            section_key = SECTION_CODES.get(code)
            if section_key:
                text_parts = []
                for elem in section.iter("text"):
                    text_parts.append("".join(elem.itertext()))
                sections[section_key] = " ".join(text_parts)[:800]
    except Exception as e:
        logger.debug("SPL XML parse error: %s", e)

    boxed_warning = sections.get("boxed_warning", "") + " " + sections.get("warnings", "")
    contraindications = sections.get("contraindications", "")
    pediatric_use = sections.get("pediatric_use", "")

    return {
        "boxed_warning": boxed_warning.strip()[:500],
        "contraindications": contraindications.strip()[:500],
        "pediatric_use": pediatric_use.strip()[:1000],
        "has_pediatric_bb_warning": _has_pediatric_bb_content(boxed_warning),
        "has_pediatric_contraindication": _has_pediatric_contraindication(contraindications + " " + pediatric_use),
        "min_age_years": _extract_min_age(pediatric_use),
        "pediatric_approved": _is_pediatric_approved(pediatric_use),
    }


def _has_pediatric_bb_content(text: str) -> bool:
    if not text:
        return False
    text_lower = text.lower()
    bb_keywords = ["pediatric", "paediatric", "children", "infant", "neonates"]
    return any(kw in text_lower for kw in bb_keywords)


def _has_pediatric_contraindication(text: str) -> bool:
    if not text:
        return False
    text_lower = text.lower()
    for pattern in _PEDIATRIC_CONTRA_PATTERNS:
        if re.search(pattern, text_lower):
            return True
    return False


def _extract_min_age(pediatric_text: str) -> float | None:
    if not pediatric_text:
        return None
    text = pediatric_text.lower()

    patterns = [
        (r"(\d+)\s+months?\s+(?:of\s+age|and\s+older)", lambda m: int(m.group(1)) / 12),
        (r"patients?\s+(\d+)\s+(?:to|through)\s+\d+\s+years?", lambda m: float(m.group(1))),
        (r"(?:age|aged)\s+(\d+)\s+years?\s+and\s+older", lambda m: float(m.group(1))),
        (r"(\d+)\s+years?\s+(?:of\s+age\s+and\s+older|and\s+older)", lambda m: float(m.group(1))),
        (r"aged?\s+(\d+)\s+(?:to\s+\d+\s+)?years?", lambda m: float(m.group(1))),
        (r"(?:age|aged)\s+(\d+)\s+year\b", lambda m: float(m.group(1))),
        (r"at\s+least\s+(\d+)\s+years?", lambda m: float(m.group(1))),
        (r"(?:neonates?|newborns?|birth|neonatal)", lambda m: 0.0),
        (r"infants?\s+(?:from|aged?)?\s*(\d+)\s+months?", lambda m: float(m.group(1)) / 12),
    ]

    ages = []
    for pattern, extractor in patterns:
        for m in re.finditer(pattern, text):
            try:
                age = extractor(m)
                ages.append(age)
            except Exception:
                pass

    return min(ages) if ages else None


def _is_pediatric_approved(pediatric_text: str) -> bool:
    if not pediatric_text:
        return False
    text = pediatric_text.lower()
    positive = [
        "is approved", "has been approved", "approved for pediatric",
        "safety and efficacy have been established",
        "safety and effectiveness have been established",
        "have been established in pediatric",
        "has been established in pediatric",
        "studied in pediatric patients",
        "is indicated for pediatric",
        "for use in pediatric patients",
        "year and older",
        "years and older",
        "years of age",
        "month of age",
        "months of age",
    ]
    negative = [
        "have not been established",
        "has not been established",
        "is not established",
        "not approved",
        "not recommended for use in patients under",
        "no data available",
        "insufficient data",
        "not studied in pediatric",
    ]
    has_positive = any(s in text for s in positive)
    has_negative = any(s in text for s in negative)
    return has_positive and not has_negative


def _empty_label() -> dict:
    return {
        "boxed_warning": "",
        "contraindications": "",
        "pediatric_use": "",
        "has_pediatric_bb_warning": False,
        "has_pediatric_contraindication": False,
        "min_age_years": None,
        "pediatric_approved": False,
    }
