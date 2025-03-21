import warnings

warnings.warn(
    "The 'SynComp' subpackage is deprecated and will be removed in future releases. "
    "Please migrate to the 'synkit' package as soon as possible,"
    + " which offers enhanced functionality. "
    "You can install it directly using pip: `pip install synkit`."
    "Code snippet: from synkit.Rule.rule_compose import RuleCompose",
    FutureWarning,
)
