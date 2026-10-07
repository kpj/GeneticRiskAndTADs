"""Modular TAD callers package and configuration grammar parser."""

import re
from typing import Any, Dict, List, Union


def expand_caller_token(token: Union[str, int]) -> List[str]:
    """Expand a caller configuration token, resolving range syntax like 'topdom_w[20..40]'.

    Syntax examples:
      - 'topdom_w[20..40]'      -> ['topdom_w20', 'topdom_w21', ..., 'topdom_w40']
      - 'topdom_w[20..40:5]'    -> ['topdom_w20', 'topdom_w25', 'topdom_w30', 'topdom_w35', 'topdom_w40']
      - 'spectraltad_l1'        -> ['spectraltad_l1']
      - 'cooltools_w250kb'      -> ['cooltools_w250kb']
    """
    token_str = str(token).strip()
    range_match = re.match(
        r"^(?P<prefix>.*)\[(?P<start>\d+)\.\.(?P<end>\d+)(?::(?P<step>\d+))?\](?P<suffix>.*)$",
        token_str,
    )
    if not range_match:
        return [token_str]

    prefix = range_match.group("prefix")
    suffix = range_match.group("suffix")
    start = int(range_match.group("start"))
    end = int(range_match.group("end"))
    step = int(range_match.group("step") or 1)

    if start > end and step > 0:
        raise ValueError(
            f"Invalid range in token '{token_str}': start ({start}) cannot exceed end ({end}) with positive step ({step})."
        )
    if step <= 0:
        raise ValueError(f"Invalid range step in token '{token_str}': step must be positive.")

    return [f"{prefix}{val}{suffix}" for val in range(start, end + 1, step)]


def expand_caller_configs(configs: Union[List[Union[str, int]], str, int]) -> List[str]:
    """Expand a list of caller configuration tokens and validate each token."""
    if isinstance(configs, (str, int)):
        configs = [configs]

    expanded: List[str] = []
    for item in configs:
        for single_token in expand_caller_token(item):
            # Validate format immediately
            parse_caller_config(single_token)
            if single_token not in expanded:
                expanded.append(single_token)
    return expanded


def resolve_consensus_caller_configs(
    selector: Union[str, List[Union[str, int]], None],
    actual_caller_configs: List[str],
) -> List[str]:
    """Resolve a consensus caller selector against a list of actual caller configurations.

    Supported selectors:
      - "all" or "*": all actual caller configs.
      - Caller tool name (e.g. "topdom", "cooltools", "spectraltad"): all actual configs matching tool.
      - Explicit config/range tokens (e.g. "topdom_w[9..11]", "cooltools_w3mb"): expanded and filtered.
      - List containing any combination of the above (e.g. ["topdom", "cooltools_w3mb"]).
    """
    if not selector:
        return []

    actual_set = set(actual_caller_configs)
    known_tools = {"topdom", "cooltools", "spectraltad", "consensustad"}

    def resolve_single_token(item: str) -> List[str]:
        item_str = str(item).strip()
        if item_str in ("all", "*"):
            return list(actual_caller_configs)

        if item_str.lower() in known_tools:
            tool_prefix = f"{item_str.lower()}_"
            return [
                c
                for c in actual_caller_configs
                if c.lower().startswith(tool_prefix) or c.lower() == item_str.lower()
            ]

        expanded = expand_caller_token(item_str)
        return [c for c in expanded if c in actual_set]

    if isinstance(selector, (str, int)):
        selectors = [selector]
    else:
        selectors = list(selector)

    result: List[str] = []
    for s in selectors:
        for matched in resolve_single_token(str(s)):
            if matched not in result:
                result.append(matched)

    return result


def parse_caller_config(caller_config: str) -> Dict[str, Any]:
    """Parse and strictly validate a caller configuration identifier into a parameter dict.

    Grammar:
      - Consensus tokens:
          'majority_vote'        -> {'caller': 'majority_vote'}
          'majority_vote_narrow' -> {'caller': 'majority_vote_narrow'}
          '*_majority_vote'      -> {'caller': '<name>'}
          '*_vote'               -> {'caller': '<name>'}
          'consensustad'         -> {'caller': 'consensustad'}

      - TopDom:
          'topdom_w<int>'        -> {'caller': 'topdom', 'window_size': int}

      - SpectralTAD:
          'spectraltad_l<int>'   -> {'caller': 'spectraltad', 'levels': int}

      - cooltools:
          'cooltools_w<int>[kb|mb]' or 'cooltools_<int>kb' -> {'caller': 'cooltools', 'window_bp': int}
    """
    token = str(caller_config).strip()

    if (
        token.endswith("_vote")
        or "vote" in token
        or token in ("majority_vote", "majority_vote_narrow", "consensustad")
    ):
        return {"caller": token}

    sep = "__" if "__" in token else "_"
    if sep not in token:
        raise ValueError(
            f"Invalid caller configuration '{token}'. Expected '<caller>_<parameter>' "
            f"(e.g. 'topdom_w20', 'spectraltad_l1', 'cooltools_w250kb', 'consensustad')."
        )

    caller, param = token.split(sep, 1)
    caller = caller.lower()

    if caller == "topdom":
        m = re.match(r"^w(?P<val>\d+)$", param)
        if not m:
            prefix_m = re.match(r"^([a-zA-Z]+)", param)
            found_prefix = prefix_m.group(1) if prefix_m else param
            raise ValueError(
                f"Invalid caller configuration '{token}': Caller 'topdom' requires parameter prefix 'w' "
                f"for window size (e.g. 'topdom_w20'), but found '{found_prefix}'."
            )
        return {"caller": "topdom", "window_size": int(m.group("val"))}

    elif caller == "spectraltad":
        m = re.match(r"^l(?P<val>\d+)$", param)
        if not m:
            prefix_m = re.match(r"^([a-zA-Z]+)", param)
            found_prefix = prefix_m.group(1) if prefix_m else param
            raise ValueError(
                f"Invalid caller configuration '{token}': Caller 'spectraltad' requires parameter prefix 'l' "
                f"for hierarchical level (e.g. 'spectraltad_l1'), but found '{found_prefix}'."
            )
        return {"caller": "spectraltad", "levels": int(m.group("val"))}

    elif caller == "cooltools":
        m = re.match(r"^(?:w)?(?P<val>\d+)(?P<unit>kb|mb|bp)?$", param.lower())
        if not m:
            raise ValueError(
                f"Invalid caller configuration '{token}': Caller 'cooltools' requires insulation window span "
                f"like 'w250kb', '250kb', or 'w250000', but found '{param}'."
            )
        val = int(m.group("val"))
        unit = m.group("unit")
        if unit == "kb":
            window_bp = val * 1_000
        elif unit == "mb":
            window_bp = val * 1_000_000
        else:
            window_bp = val
        return {"caller": "cooltools", "window_bp": window_bp}

    else:
        raise ValueError(
            f"Unknown TAD caller '{caller}' in configuration '{token}'. "
            f"Supported callers: 'topdom', 'spectraltad', 'cooltools', 'consensustad'."
        )
