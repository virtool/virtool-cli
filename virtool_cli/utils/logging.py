import os
import sys
import structlog
import logging

logging.basicConfig(
    format="%(message)s",
    stream=sys.stderr,
    level=logging.DEBUG,
)

shared_processors = [
    structlog.stdlib.filter_by_level,
    structlog.stdlib.add_logger_name,
    structlog.stdlib.add_log_level,
    structlog.stdlib.PositionalArgumentsFormatter(),
    structlog.processors.StackInfoRenderer(),
    structlog.processors.UnicodeDecoder(),
    structlog.processors.CallsiteParameterAdder(
            {
                structlog.processors.CallsiteParameter.FILENAME,
                structlog.processors.CallsiteParameter.FUNC_NAME,
                structlog.processors.CallsiteParameter.LINENO,
            }
        ),
]
if sys.stderr.isatty() and not os.environ.get("NO_COLOR"):
    processors = shared_processors + [
        structlog.dev.ConsoleRenderer()
    ]
else:
    processors = shared_processors + [
        structlog.processors.TimeStamper(fmt="iso"),
        structlog.processors.format_exc_info,
        structlog.processors.dict_tracebacks,
        structlog.processors.JSONRenderer(),
    ]

structlog.configure(
    processors=processors,
    logger_factory=structlog.stdlib.LoggerFactory(),
    cache_logger_on_first_use=True,
)

base_logger = structlog.get_logger()

DEFAULT_LOGGER = structlog.make_filtering_bound_logger(logging.INFO)
DEBUG_LOGGER = structlog.make_filtering_bound_logger(logging.DEBUG)
