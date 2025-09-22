Hello! I need your help debugging some failing tests in my Julia project.

Please follow these steps:

1.  **Run the tests**: Execute the specific test file using the following command.
    ```sh
    julia --project=. -e 'include("${TEST_FILE}")'
    ```

2.  **Identify Failures**: Review the output and identify all the tests that are failing.

3.  **Analyze and Diagnose**: For each failing test, please investigate both the test script and the underlying source code it is testing. Determine the root cause of the failure:
    *   **Is the function correct and the test is outdated or incorrect?**
    *   **Is the test correct and the function has a bug?**
    *   **Is the test correct and the function is just a placeholder that needs full implementation?**

4.  **Propose a Solution**: Based on your diagnosis, please provide the necessary code changes.
    *   If the **function is correct**, update the test to match the expected behavior.
    *   If the **function is broken or incomplete**, please fix the implementation so that the test passes.