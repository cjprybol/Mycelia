"""
Utils Submodule for Mycelia.jl

Handles utility functions, external tool integration, and system management.
"""
module Utils

import Conda
import HTTP
import JSON
import ProgressMeter
import CSV
import Dates
import Statistics
import Random
import SHA
import UUIDs

# ============================================================================
# PUBLIC API EXPORTS
# ============================================================================

# General utilities
export normalized_current_datetime, check_memory_usage
export estimate_memory_requirements, validate_file_paths

# Bioconda integration
export add_bioconda_env, setup_bioconda
export install_tool, check_tool_availability

# SLURM integration
export submit_slurm_job, check_job_status
export slurm_array_job, cancel_job

# Cloud storage (Rclone)
export setup_rclone, sync_data
export upload_results, download_datasets

# Container integration
export run_docker_container, run_singularity_container
export build_container, manage_containers

# File and system utilities
export safe_file_operations, backup_files
export cleanup_temp_files, monitor_disk_space

# ============================================================================
# INTERNAL HELPER FUNCTIONS (will get _ prefix)
# ============================================================================

export _validate_tool_installation, _setup_environment
export _parse_job_output, _handle_tool_errors
export _system_checks, _resource_monitoring

# ============================================================================
# FUNCTION IMPLEMENTATIONS
# ============================================================================

# Functions will be moved from:
# - utility-functions.jl
# - bioconda.jl
# - slurm-sbatch.jl  
# - rclone.jl
# to here

end # module Utils