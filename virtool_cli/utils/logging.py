import sys
import structlog

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
            }
        ),
]
if sys.stderr.isatty():
    processors = shared_processors + [
        structlog.dev.ConsoleRenderer(),
        structlog.processors.CallsiteParameterAdder(
            {
                structlog.processors.CallsiteParameter.LINENO,
            }
        )
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
    wrapper_class=structlog.stdlib.BoundLogger,
    logger_factory=structlog.stdlib.LoggerFactory(),
    cache_logger_on_first_use=True,
)

base_logger = structlog.get_logger()