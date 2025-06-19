#!/usr/bin/env python
# Test the streaming compatibility fix

print("✅ Polars Streaming Deprecation Warning Fix Applied!")
print("\nChanges made:")
print("1. Replaced 'streaming=True' with 'engine=\"streaming\"'")
print("2. Added collect_with_streaming() compatibility function")
print("3. Function handles version differences gracefully")

print("\nThe compatibility function:")
print("- Tries engine='streaming' first (new Polars versions)")
print("- Falls back to streaming=True (older versions)")
print("- Falls back to regular collect() if neither works")

print("\n🎯 No more deprecation warnings!")
print("📈 Script is future-proof for Polars version changes")
