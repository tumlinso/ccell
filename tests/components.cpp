#include <cmath>
#include <type_traits>

import ccell.anndata;
import ccell.preprocess;
import ccell.velocity;
import anndatapp;
import preprocess;
import velocitypp;
import velocitypp.precision;

static_assert(std::is_same_v<ccell::anndata::Array, anndatapp::Array>);
static_assert(std::is_same_v<ccell::preprocess::QcStatus, preprocess::QcStatus>);
static_assert(std::is_same_v<ccell::velocity::PreparedInput, velocitypp::PreparedInput>);

#include "workflow.hh"

int main() { return workflow(); }
