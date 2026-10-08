// Record of changes to the calculation engines that move results. Every entry says what changed, why,
// how much typical results moved, and how to reproduce the earlier behaviour where that is possible.
// Results, reports and case files are stamped with the build that produced them, so a number can always be traced.
// Entry shape: { date: 'yyyy-mm-dd', suites: ['flow', …], title, what, effect, revert }.
export const CHANGELOG = [];
/** Entries that concern one suite, newest first. */
export const changesFor = (id) => CHANGELOG.filter((c) => c.suites.includes(id));
