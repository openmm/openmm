# AI Policy

This document describes the policy on using artificial intelligence in OpenMM development.  This field is changing very
quickly, and this policy is expected to evolve with time.

## What Is Required

Any use of AI in issues and pull requests must be disclosed.  Examples might include using AI to fix a bug, implement
a new feature, or generate an issue summary.  In every case it must be stated that AI was used.  This requirement is not
meant to discourage any particular use of AI.  It is simply so that everyone understands what they are dealing with.

## What AI Is Good For

### Finding Bugs

AI can be a powerful tool for identifying bugs and understanding their causes.  We welcome all bug reports, whatever
tools you used in investigating them.

### Fixing Bugs

In some cases, AI can also be useful for suggesting how to fix bugs.  This is mainly true for localized fixes that
change only a few lines of code and are easy to verify.  It is generally less suitable for more involved fixes that
affect many pieces of code or involve large amounts of new code.

### Writing Code That Will Not Be Included In OpenMM

For example, when reporting a bug, it is always best if you can provide a self-contained test that reproduces the
problem.  You should feel free to use AI in creating the test, as long as you verify that it runs and really does
reproduce the bug.

## What AI Is Not Good For

### Writing New Code

We do not have a strict policy against AI generated code, but given the current state of the field, it usually is not
suitable for this purpose.

AI generated code is held to the same standards as any other submitted code.  It should be clean, correct, efficient,
easy to understand, and easy to maintain.  In our experience, AI generated code tends to be the opposite: convoluted,
inefficient, hard to understand, inconsistent with the rest of the code base, based on incorrect assumptions, and filled
with subtly incorrect behaviors.

You are responsible for the correctness of any code you submit.  For AI generated code, that means you should have
carefully reviewed it to ensure it is correct.  You also are expected to fully understand how it works and be prepared
to answer questions about it.

You also must confirm that you have the legal right to submit it.  Most AI coding models are trained on large amounts of
copyrighted code used without permission.  They are able to exactly or almost exactly reproduce large blocks of the code
they were trained on.  Laundering copyrighted code through an AI model does not remove the copyright or free you from
the obligation to follow the terms of the license.  Any large block of AI generated code therefore has a risk of being
encumbered by copyright.  If you are not absolutely certain you have the legal right to submit it, please do not submit
it!

### Replying To Questions

When someone asks a question about an issue or pull request, they are asking you, not an AI.  Please do not use AI to
generate an answer!